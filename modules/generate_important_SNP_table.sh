#!/bin/bash
#SBATCH --job-name=ALS_variant_evidence
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/generate_important_SNP_table.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/generate_important_SNP_table.err
#SBATCH --time=06:00:00
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail

ROOT="$HOME/donglab/data/target_ALS"
GWAS="$HOME/donglab/references/GWAS/ALS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
GTF="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf"
OUTDIR="$ROOT/QTL/ALS_variant_validation_evidence"

mkdir -p "$OUTDIR"

for f in "$GWAS" "$GTF"; do
    [[ -s "$f" ]] || { echo "ERROR: missing/empty $f"; exit 1; }
done

python -u - "$GWAS" "$ROOT" "$OUTDIR" "$GTF" <<'PY'
import sys
import csv
import os
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

GWAS=Path(sys.argv[1])
ROOT=Path(sys.argv[2])
OUTDIR=Path(sys.argv[3])
GTF=Path(sys.argv[4])

TISSUES=[
    "Cortex",
    "Spinal_Cord",
    "Thoracic_Spinal_Cord",
    "Cerebellum",
    "Frontal_Cortex",
    "Cervical_Spinal_Cord",
    "Lumbar_Spinal_Cord",
    "Motor_Cortex"
]

QTLS=["eQTL","sQTL","cQTL"]

ALS_GENES = {x.upper() for x in [
    'ATXN1', 'ATXN2', 'HOMER1', 'C9ORF72', 'SOD1', 'FUS',
    'STMN2', 'TARDBP', 'TBK1', 'UNC13A', 'ACSL5', 'ALS2',
    'ANG', 'ANXA11', 'C21ORF2', 'CAV1', 'CCNF', 'CHCHD10',
    'CHMP2B', 'CHRNA3', 'DAO', 'DCTN1', 'DNAJC7', 'ELP3',
    'ERBB4', 'EWSR1', 'FIG4', 'GLE1', 'GLT8D1', 'HNRNPA1',
    'HNRNPA2B1', 'KANK1', 'KIF5A', 'LGALSL', 'MATR3', 'MOBP',
    'NEFH', 'NEK1', 'NIPA1', 'OPTN', 'PARK9', 'PFN1',
    'PON1', 'PON2', 'PON3', 'PRPH', 'SARM1', 'SCFD1',
    'SETX', 'SIGMAR1', 'SPG11', 'SPTLC1', 'SQSTM1', 'TAF15',
    'TIA1', 'TUBA4A', 'UBQLN2', 'VAPB', 'VCP', 'WDR7',
]}

NEUROLOGY_GENES = {x.upper() for x in [
    'RIMS1', 'RIMS2', 'SNCA', 'SMN1', 'FMR1', 'PRKN', 'PINK1',
    'PARK7', 'APP', 'PSEN1', 'PSEN2', 'APOE', 'LRRK2', 'PICALM',
    'ABCA7', 'CLU', 'CR1', 'PLD3', 'TREM2', 'SORL1', 'PARK2',
    'GCH1', 'TH', 'ATP13A2', 'VPS35', 'HTT', 'ALSIN', 'TUBB4A', 'TOR1A', 'FXN', 'SIL1',
]}
# Keep user labels; resolve aliases through the same GENCODE gene IDs where possible.
RBP_FILE = Path.home() / "donglab/references/genome/Homo_sapiens/RBP2GO/RNA_binding_proteins/RNABP_list.csv"
if not RBP_FILE.is_file():
    raise FileNotFoundError(f"Missing RNA-binding protein list: {RBP_FILE}")
rbp_table = pd.read_csv(RBP_FILE, comment="#", dtype=str)
if "Gene_Name" not in rbp_table.columns:
    raise RuntimeError("RNA-binding protein list requires a Gene_Name column")
RBP_GENES = set(rbp_table["Gene_Name"].dropna().str.strip().str.upper())
RBP_GENES.discard("")

# Gene-list points use only targets in the significant QTL input files.
# Keep these separate from candidate genes pooled across SMR/coloc results.
significant_qtl_gene_symbols = defaultdict(set)


# ------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------

def norm(x):
    return re.sub(r"[^a-z0-9]","",str(x).lower())

def strip_gene_version(x):
    return re.sub(r"\.\d+$","",str(x))

def chr_key(x):
    return re.sub(r"^chr","",str(x),flags=re.I).upper()

def clean_value(x):
    if x is None:
        return None
    x=str(x).strip()
    if not x or x.lower() in {"nan","na","none",".","null"}:
        return None
    return x

def clean_snp(x):
    return clean_value(x)

def find_col(cols,aliases,required=False):
    lookup={norm(c):c for c in cols}

    for a in aliases:
        if norm(a) in lookup:
            return lookup[norm(a)]

    for c in cols:
        nc=norm(c)
        for a in aliases:
            na=norm(a)
            if na and na in nc:
                return c

    if required:
        raise RuntimeError(
            f"Could not find any of {aliases}; columns={list(cols)}"
        )

    return None

def read_table(path):
    try:
        return pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            low_memory=False
        )
    except Exception:
        return pd.read_csv(
            path,
            sep=r"\s+",
            dtype=str,
            engine="python"
        )

def parse_gtf_attrs(s):
    out={}
    for item in str(s).strip().split(";"):
        item=item.strip()
        if not item:
            continue

        parts=item.split(" ",1)
        if len(parts)==2:
            out[parts[0]]=parts[1].strip().strip('"')

    return out

def parse_event_coordinates(event):
    """
    Supported event formats:

    sQTL:
      chr10:+:100381448-100394498

    cQTL:
      chr10:100164003-100167406:-

    legacy:
      chr10:100164003:100167406:-

    no-strand:
      chr10:100164003-100167406
      chr10:100164003:100167406
    """
    s=str(event).strip()

    patterns=[
        (r"^(chr[^:]+):([+-]):(\d+)-(\d+)$","strand_first"),
        (r"^(chr[^:]+):(\d+)-(\d+):([+-])$","strand_last"),
        (r"^(chr[^:]+):(\d+):(\d+):([+-])$","strand_last"),
        (r"^(chr[^:]+):(\d+)-(\d+)$","no_strand"),
        (r"^(chr[^:]+):(\d+):(\d+)$","no_strand")
    ]

    for pat,mode in patterns:
        m=re.match(pat,s)
        if not m:
            continue

        if mode=="strand_first":
            chrom,strand,start,end=m.groups()

        elif mode=="strand_last":
            chrom,start,end,strand=m.groups()

        else:
            chrom,start,end=m.groups()
            strand=None

        start=int(start)
        end=int(end)

        if start>end:
            start,end=end,start

        return {
            "chr":chrom,
            "start":start,
            "end":end,
            "strand":strand
        }

    return None

# ------------------------------------------------------------------
# LOAD GENCODE GENE ANNOTATION
# ------------------------------------------------------------------

print(f"Loading GENCODE genes: {GTF}",flush=True)

GENES_BY_CHR=defaultdict(list)
GENE_BY_SYMBOL={}
GENE_BY_ID={}

with open(GTF) as f:
    for line in f:
        if line.startswith("#"):
            continue

        x=line.rstrip("\n").split("\t")

        if len(x)<9 or x[2]!="gene":
            continue

        attrs=parse_gtf_attrs(x[8])

        gene_id=attrs.get("gene_id","")
        gene_symbol=attrs.get("gene_name","")

        if not gene_id:
            continue

        g={
            "chr":x[0],
            "start":int(x[3])-1,
            "end":int(x[4]),
            "strand":x[6],
            "gene_id":gene_id,
            "gene_id_base":strip_gene_version(gene_id),
            "gene_symbol":gene_symbol
        }

        GENES_BY_CHR[chr_key(x[0])].append(g)

        GENE_BY_ID[g["gene_id_base"].upper()]=g

        if gene_symbol:
            GENE_BY_SYMBOL[norm(gene_symbol)]=g

for chrom in GENES_BY_CHR:
    GENES_BY_CHR[chrom].sort(
        key=lambda g:(g["start"],g["end"])
    )

print(
    f"Loaded {sum(len(v) for v in GENES_BY_CHR.values()):,} GENCODE genes",
    flush=True
)

def genes_for_trait(trait,qtl_type):
    """
    Return GENCODE genes associated with a QTL trait.

    eQTL:
      direct gene-symbol or Ensembl-ID lookup

    sQTL/cQTL:
      any GENCODE gene whose genomic span overlaps the event by >=1 bp

    Strand is deliberately NOT required for sQTL/cQTL overlap.
    """
    trait=clean_value(trait)

    if trait is None:
        return []

    if qtl_type=="eQTL":
        base=strip_gene_version(trait).upper()

        if base in GENE_BY_ID:
            return [GENE_BY_ID[base]]

        n=norm(trait)

        if n in GENE_BY_SYMBOL:
            return [GENE_BY_SYMBOL[n]]

        return []

    loc=parse_event_coordinates(trait)

    if loc is None:
        return []

    chrom=chr_key(loc["chr"])
    start=loc["start"]
    end=loc["end"]

    hits=[]

    for g in GENES_BY_CHR.get(chrom,[]):
        if g["start"] > end:
            break

        if start <= g["end"] and end >= g["start"]:
            hits.append(g)

    return hits

# ------------------------------------------------------------------
# 1. LOAD ALS GWAS
# ------------------------------------------------------------------

print(f"Loading GWAS: {GWAS}",flush=True)

gwas=read_table(GWAS)

snp_col=find_col(
    gwas.columns,
    ["rsid","rs_id","variant_id","variant","snp","snpid"],
    required=True
)

p_col=find_col(
    gwas.columns,
    ["p_value","pvalue","pval","p"]
)

beta_col=find_col(
    gwas.columns,
    ["beta","effect","effect_size"]
)

se_col=find_col(
    gwas.columns,
    ["standard_error","stderr","se"]
)

ea_col=find_col(
    gwas.columns,
    ["effect_allele","effectallele","ea"]
)

oa_col=find_col(
    gwas.columns,
    ["other_allele","otherallele","nea","non_effect_allele"]
)

chr_col=find_col(
    gwas.columns,
    ["chromosome","chrom","chr"]
)

pos_col=find_col(
    gwas.columns,
    ["base_pair_location","position","pos","bp"]
)

gwas["SNP"]=gwas[snp_col].map(clean_snp)

gwas=(
    gwas[
        gwas["SNP"].notna()
    ]
    .drop_duplicates("SNP")
    .copy()
)

out=pd.DataFrame({
    "SNP":gwas["SNP"]
})

def copy_gwas(col,new_name):
    if col:
        out[new_name]=gwas[col].values
    else:
        out[new_name]=np.nan

copy_gwas(chr_col,"GWAS_chr")
copy_gwas(pos_col,"GWAS_pos")
copy_gwas(ea_col,"GWAS_effect_allele")
copy_gwas(oa_col,"GWAS_other_allele")
copy_gwas(beta_col,"GWAS_beta")
copy_gwas(se_col,"GWAS_se")
copy_gwas(p_col,"GWAS_p")

out["GWAS_p"]=pd.to_numeric(
    out["GWAS_p"],
    errors="coerce"
)

# USER-REQUESTED GWAS EVIDENCE:
# nominal GWAS significance p < 0.05
out["GWAS_sig"]=(
    out["GWAS_p"].notna()
    &
    (out["GWAS_p"] < 0.05)
)

out=out.set_index(
    "SNP",
    drop=False
)

print(
    f"GWAS variants loaded: {len(out):,}",
    flush=True
)

print(
    f"GWAS p<0.05: {int(out['GWAS_sig'].sum()):,}",
    flush=True
)

# ------------------------------------------------------------------
# EVIDENCE ACCUMULATORS
# ------------------------------------------------------------------

gene_symbols=defaultdict(set)
gene_ids=defaultdict(set)
gene_pairs=defaultdict(set)

def add_gene(snp,g):
    symbol=clean_value(g.get("gene_symbol"))
    gid=clean_value(g.get("gene_id"))

    if symbol:
        gene_symbols[snp].add(symbol)

    if gid:
        gene_ids[snp].add(gid)

    if symbol and gid:
        gene_pairs[snp].add(f"{symbol}|{gid}")
    elif gid:
        gene_pairs[snp].add(gid)
    elif symbol:
        gene_pairs[snp].add(symbol)

for q in QTLS:
    out[f"{q}_significant"]=False
    out[f"{q}_tissues"]=""
    out[f"{q}_traits"]=""
    out[f"{q}_best_FDR"]=np.nan

    out[f"{q}_SMR_significant"]=False
    out[f"{q}_SMR_tissues"]=""
    out[f"{q}_best_p_SMR"]=np.nan
    out[f"{q}_HEIDI_pass"]=False

    out[f"{q}_coloc_significant"]=False
    out[f"{q}_coloc_tissues"]=""
    out[f"{q}_best_PPH4"]=np.nan

    out[f"{q}_SuSiE_credible_set"]=False
    out[f"{q}_best_PIP"]=np.nan

qtl_tissues={
    q:defaultdict(set)
    for q in QTLS
}

qtl_traits={
    q:defaultdict(set)
    for q in QTLS
}

smr_tissues={
    q:defaultdict(set)
    for q in QTLS
}

coloc_tissues={
    q:defaultdict(set)
    for q in QTLS
}

# ------------------------------------------------------------------
# 2. SIGNIFICANT QTL RESULTS
# ------------------------------------------------------------------

for tissue in TISSUES:
    for q in QTLS:

        f=(
            ROOT
            / tissue
            / q
            / "results"
            / f"{tissue}_{q}.FDR0.05.txt"
        )

        if not f.exists():
            print(
                f"MISSING QTL: {f}",
                flush=True
            )
            continue

        print(
            f"Reading {q}: {tissue}",
            flush=True
        )

        d=read_table(f)

        sc=find_col(
            d.columns,
            ["snpid","snp_id","snp","rsid","variant"],
            required=True
        )

        fc=find_col(
            d.columns,
            ["fdr","qvalue","q_value","padj"]
        )

        tc=find_col(
            d.columns,
            [
                "geneid",
                "gene_id",
                "junction_id",
                "junctionid",
                "circ_id",
                "circid",
                "phenotype_id",
                "probeid",
                "event_id"
            ]
        )

        gc=find_col(
            d.columns,
            ["gene_symbol","gene_name","symbol"]
        )

        for _,r in d.iterrows():

            snp=clean_snp(
                r.get(sc)
            )

            if snp not in out.index:
                continue

            out.at[
                snp,
                f"{q}_significant"
            ]=True

            qtl_tissues[q][snp].add(
                tissue
            )

            trait=None

            if tc:
                trait=clean_value(
                    r.get(tc)
                )

            if trait:
                qtl_traits[q][snp].add(
                    trait
                )

                # Map every QTL trait back to GENCODE.
                for g in genes_for_trait(
                    trait,
                    q
                ):
                    add_gene(
                        snp,
                        g
                    )
                    if g.get("gene_symbol"):
                        significant_qtl_gene_symbols[snp].add(g["gene_symbol"])

            # Preserve any explicit gene-symbol column too.
            if gc:
                symbol=clean_value(
                    r.get(gc)
                )

                if symbol:
                    significant_qtl_gene_symbols[snp].add(symbol)
                    g=GENE_BY_SYMBOL.get(
                        norm(symbol)
                    )

                    if g:
                        add_gene(
                            snp,
                            g
                        )
                    else:
                        gene_symbols[snp].add(
                            symbol
                        )

            if fc:
                try:
                    x=float(
                        r[fc]
                    )

                    old=out.at[
                        snp,
                        f"{q}_best_FDR"
                    ]

                    if (
                        pd.isna(old)
                        or x<old
                    ):
                        out.at[
                            snp,
                            f"{q}_best_FDR"
                        ]=x

                except Exception:
                    pass

# ------------------------------------------------------------------
# 3. SMR + HEIDI
# ------------------------------------------------------------------

for tissue in TISSUES:
    for q in QTLS:

        f=(
            ROOT
            / tissue
            / "MR"
            / f"{q}_SMR_HEIDI"
            / "results"
            / f"{tissue}_{q}_SMR_HEIDI.smr"
        )

        if not f.exists():
            print(
                f"MISSING SMR: {f}",
                flush=True
            )
            continue

        print(
            f"Reading SMR: {q} {tissue}",
            flush=True
        )

        d=read_table(f)

        sc=find_col(
            d.columns,
            [
                "topSNP",
                "top_snp",
                "snp",
                "rsid"
            ]
        )

        pc=find_col(
            d.columns,
            [
                "p_SMR",
                "psmr"
            ]
        )

        hc=find_col(
            d.columns,
            [
                "p_HEIDI",
                "pheidi"
            ]
        )

        gc=find_col(
            d.columns,
            [
                "Gene",
                "gene",
                "gene_symbol"
            ]
        )

        probe_col=find_col(
            d.columns,
            [
                "probeID",
                "probe_id",
                "phenotype_id",
                "geneid",
                "gene_id",
                "junction_id",
                "circ_id"
            ]
        )

        if sc is None:
            continue

        # ----------------------------------------------------------
        # Recalculate tissue-level BH FDR from CURRENT .smr results
        # ----------------------------------------------------------

        fdr=pd.Series(
            np.nan,
            index=d.index
        )

        if pc:
            p=pd.to_numeric(
                d[pc],
                errors="coerce"
            )

            valid=p.notna()

            if valid.any():
                vals=p[valid].sort_values()

                m=len(vals)

                adj=(
                    vals
                    * m
                    / np.arange(
                        1,
                        m+1
                    )
                )

                adj=adj.clip(
                    upper=1
                )

                adj=np.minimum.accumulate(
                    adj.iloc[::-1].values
                )[::-1]

                fdr.loc[
                    vals.index
                ]=adj

        for idx,r in d.iterrows():

            snp=clean_snp(
                r.get(sc)
            )

            if snp not in out.index:
                continue

            # Map SMR gene/probe to GENCODE where possible.
            smr_trait=None

            if probe_col:
                smr_trait=clean_value(
                    r.get(probe_col)
                )

            if smr_trait:
                for g in genes_for_trait(
                    smr_trait,
                    q
                ):
                    add_gene(
                        snp,
                        g
                    )

            if gc:
                gene_label=clean_value(
                    r.get(gc)
                )

                if gene_label:
                    for g in genes_for_trait(
                        gene_label,
                        q
                    ):
                        add_gene(
                            snp,
                            g
                        )

            # Best raw SMR p-value.
            try:
                psmr=(
                    float(r[pc])
                    if pc
                    else np.nan
                )

                old=out.at[
                    snp,
                    f"{q}_best_p_SMR"
                ]

                if (
                    np.isfinite(psmr)
                    and (
                        pd.isna(old)
                        or psmr<old
                    )
                ):
                    out.at[
                        snp,
                        f"{q}_best_p_SMR"
                    ]=psmr

            except Exception:
                pass

            # Tissue-level SMR FDR < 0.05.
            try:
                if float(
                    fdr.loc[idx]
                ) < 0.05:

                    out.at[
                        snp,
                        f"{q}_SMR_significant"
                    ]=True

                    smr_tissues[q][snp].add(
                        tissue
                    )

            except Exception:
                pass

            # HEIDI pass = p_HEIDI >= 0.05.
            if hc:
                try:
                    if float(
                        r[hc]
                    ) >= 0.05:

                        out.at[
                            snp,
                            f"{q}_HEIDI_pass"
                        ]=True

                except Exception:
                    pass

# ------------------------------------------------------------------
# 4. SuSiE / COLOC
#
# Search recursively because the three pipelines write multiple per-locus/global TSV files under each *_SuSiE_coloc directory.
# ------------------------------------------------------------------

for tissue in TISSUES:
    for q in QTLS:

        droot=(
            ROOT
            / tissue
            / "MR"
            / f"{q}_SuSiE_coloc"
        )

        if not droot.exists():
            continue

        print(
            f"Scanning SuSiE/coloc: {q} {tissue}",
            flush=True
        )

        for f in droot.rglob("*.tsv*"):

            try:
                d=read_table(f)
            except Exception:
                continue

            sc=find_col(
                d.columns,
                [
                    "candidate_snp",
                    "snp",
                    "snpid",
                    "rsid",
                    "variant"
                ]
            )

            if sc is None:
                continue

            pph4c=find_col(
                d.columns,
                [
                    "PPH4",
                    "PP.H4",
                    "H4"
                ]
            )

            h4pass=find_col(
                d.columns,
                [
                    "H4_reference_pass",
                    "H4_pass",
                    "coloc_pass"
                ]
            )

            pipc=find_col(
                d.columns,
                [
                    "PIP",
                    "susie_pip",
                    "posterior_inclusion_probability"
                ]
            )

            csc=find_col(
                d.columns,
                [
                    "credible_set",
                    "in_credible_set",
                    "cs"
                ]
            )

            trait_col=find_col(
                d.columns,
                [
                    "geneid",
                    "gene_id",
                    "gene",
                    "gene_symbol",
                    "probeID",
                    "probe_id",
                    "junction_id",
                    "junctionid",
                    "circ_id",
                    "circid",
                    "phenotype_id",
                    "event_id"
                ]
            )

            # Ignore unrelated SNP-containing tables.
            if not any([
                pph4c,
                h4pass,
                pipc,
                csc
            ]):
                continue

            for _,r in d.iterrows():

                snp=clean_snp(
                    r.get(sc)
                )

                if snp not in out.index:
                    continue

                # Attach candidate gene annotation from coloc trait.
                if trait_col:
                    trait=clean_value(
                        r.get(trait_col)
                    )

                    if trait:
                        for g in genes_for_trait(
                            trait,
                            q
                        ):
                            add_gene(
                                snp,
                                g
                            )

                # Best PPH4.
                if pph4c:
                    try:
                        x=float(
                            r[pph4c]
                        )

                        old=out.at[
                            snp,
                            f"{q}_best_PPH4"
                        ]

                        if (
                            pd.isna(old)
                            or x>old
                        ):
                            out.at[
                                snp,
                                f"{q}_best_PPH4"
                            ]=x

                    except Exception:
                        pass

                # Prefer explicit H4_reference_pass.
                passed=False

                if h4pass:
                    v=str(
                        r[h4pass]
                    ).strip().lower()

                    passed=(
                        v
                        in {
                            "true",
                            "t",
                            "1",
                            "yes",
                            "pass"
                        }
                    )

                # Fallback only if explicit pass column absent.
                elif pph4c:
                    try:
                        passed=(
                            float(
                                r[pph4c]
                            )
                            >=0.80
                        )
                    except Exception:
                        pass

                if passed:
                    out.at[
                        snp,
                        f"{q}_coloc_significant"
                    ]=True

                    coloc_tissues[q][snp].add(
                        tissue
                    )

                # Best SuSiE PIP.
                if pipc:
                    try:
                        x=float(
                            r[pipc]
                        )

                        old=out.at[
                            snp,
                            f"{q}_best_PIP"
                        ]

                        if (
                            pd.isna(old)
                            or x>old
                        ):
                            out.at[
                                snp,
                                f"{q}_best_PIP"
                            ]=x

                    except Exception:
                        pass

                # Credible-set membership.
                if csc:
                    v=str(
                        r[csc]
                    ).strip().lower()

                    is_member = v in {"true", "t", "yes", "y"}
                    if not is_member:
                        try:
                            cs_number = float(v)
                            is_member = np.isfinite(cs_number) and cs_number > 0 and cs_number.is_integer()
                        except (TypeError, ValueError):
                            pass
                    if is_member:
                        out.at[
                            snp,
                            f"{q}_SuSiE_credible_set"
                        ]=True

# ------------------------------------------------------------------
# 5. COLLAPSE SETS
# ------------------------------------------------------------------

for q in QTLS:

    for snp,v in qtl_tissues[q].items():
        out.at[
            snp,
            f"{q}_tissues"
        ]=";".join(
            sorted(v)
        )

    for snp,v in qtl_traits[q].items():
        out.at[
            snp,
            f"{q}_traits"
        ]=";".join(
            sorted(v)
        )

    for snp,v in smr_tissues[q].items():
        out.at[
            snp,
            f"{q}_SMR_tissues"
        ]=";".join(
            sorted(v)
        )

    for snp,v in coloc_tissues[q].items():
        out.at[
            snp,
            f"{q}_coloc_tissues"
        ]=";".join(
            sorted(v)
        )

# GENCODE candidate-gene annotations.
out["candidate_gene_symbols"]=[
    ";".join(
        sorted(
            gene_symbols.get(
                snp,
                set()
            )
        )
    )
    for snp in out.index
]

out["candidate_gene_ids"]=[
    ";".join(
        sorted(
            gene_ids.get(
                snp,
                set()
            )
        )
    )
    for snp in out.index
]

out["candidate_genes"]=[
    ";".join(
        sorted(
            gene_pairs.get(
                snp,
                set()
            )
        )
    )
    for snp in out.index
]

# Include current symbols for user-provided alternate gene labels.
# GENCODE may use ATP13A2 for PARK9 and CFAP410 for C21orf2.
ALS_GENE_MATCH_NAMES = ALS_GENES | {"ATP13A2", "CFAP410"}
# Gene-list matches are case insensitive; retain original symbols in output.
# Include known alternate labels for the same neurology-panel genes.
NEUROLOGY_GENE_MATCH_NAMES = NEUROLOGY_GENES | {"PARK9"}
for label, gene_list in [("ALS_relevant_gene", ALS_GENE_MATCH_NAMES),
                         ("NEUROLOGY_relevant_gene", NEUROLOGY_GENE_MATCH_NAMES),
                         ("RNA_binding_protein", RBP_GENES)]:
    matches = [
        sorted(s for s in significant_qtl_gene_symbols.get(snp, set())
               if s.upper() in gene_list)
        for snp in out.index
    ]
    out[f"{label}_matches"] = [";".join(v) for v in matches]
    out[label] = [bool(v) for v in matches]

# Observed homozygote differences, using the supplied ranking-script rules.
# Optional TSV named by QTL_MATRIX_CONFIG: tissue, qtl_type, matrix_path, scale.
# scale must be zscore (eQTL), fraction or percent (sQTL/cQTL).
# Paths may contain {tissue} and {root}. sQTL has a confirmed default.
MIN_HOMOZYGOTES = 5
THRESHOLDS = {"eQTL": 0.5, "sQTL": 0.5, "cQTL": 0.1}
config = {}
config_file = os.environ.get("QTL_MATRIX_CONFIG")
if config_file:
    cfg = pd.read_csv(Path(config_file).expanduser(), sep="\t", dtype=str, keep_default_na=False)
    required = {"tissue", "qtl_type", "matrix_path", "scale"}
    if not required.issubset(cfg.columns):
        raise ValueError(f"Matrix configuration requires columns: {sorted(required)}")
    for r in cfg.to_dict("records"):
        key = (r["tissue"], r["qtl_type"])
        if key in config:
            raise ValueError(f"Duplicate matrix configuration: {key}")
        if key[0] not in TISSUES or key[1] not in QTLS:
            raise ValueError(f"Unknown tissue/QTL in configuration: {key}")
        allowed = {"zscore"} if key[1] == "eQTL" else {"fraction", "percent"}
        if r["scale"] not in allowed:
            raise ValueError(f"Invalid scale for {key}: expected {allowed}")
        config[key] = (Path(r["matrix_path"].format(root=ROOT, tissue=key[0])).expanduser(), r["scale"])

def selected_matrix(path, wanted):
    values = {}
    with path.open() as handle:
        samples = handle.readline().rstrip("\r\n").split("\t")[1:]
        if not samples or len(samples) != len(set(samples)):
            raise ValueError(f"Empty/duplicate sample IDs: {path}")
        for line in handle:
            key, _, tail = line.rstrip("\r\n").partition("\t")
            if key not in wanted:
                continue
            if key in values:
                raise ValueError(f"Duplicate matrix row {key}: {path}")
            v = pd.to_numeric(pd.Series(tail.split("\t")), errors="coerce").to_numpy(dtype=float)
            if len(v) != len(samples):
                raise ValueError(f"Wrong matrix row length {key}: {path}")
            values[key] = v
    return samples, values

jobs = []
audit = []
def audit_failure(t, q, pairs, reason):
    for snp, trait in pairs:
        audit.append(dict(SNP=snp, tissue=t, qtl_type=q, trait=trait, status=reason, passes=False))

# Use all significant pairs, not only the top-for-boxplot subset.
for t in TISSUES:
    for q in QTLS:
        # Re-read significant rows to retain exact tissue/SNP/trait pairs.
        sigpath = ROOT/t/q/"results"/f"{t}_{q}.FDR0.05.txt"
        if not sigpath.exists():
            continue
        sig = read_table(sigpath)
        sc = find_col(sig.columns, ["snpid", "snp_id", "snp", "rsid", "variant"], required=True)
        tc = find_col(sig.columns, ["geneid", "gene_id", "junction_id", "junctionid", "circ_id", "circid", "phenotype_id", "probeid", "event_id"])
        if tc is None:
            pairs = sorted({(clean_snp(v), "") for v in sig[sc] if clean_snp(v) in out.index})
            audit_failure(t, q, pairs, "missing_trait_column")
            continue
        pairs = sorted({(clean_snp(r[sc]), clean_value(r[tc])) for r in sig.to_dict("records")
                        if clean_snp(r[sc]) in out.index and clean_value(r[tc]) is not None})
        if not pairs:
            continue
        setting = config.get((t, q))
        if setting is None and q == "sQTL":
            setting = (ROOT/t/q/f"splicing_{t}.txt", "fraction")
        if setting is None:
            audit_failure(t, q, pairs, "matrix_configuration_required")
            continue
        matrix_path, scale = setting
        genotype_path = ROOT/t/q/f"snp_{t}.txt"
        locpath = ROOT/t/q/"snp_location.txt"
        if not all(f.is_file() for f in (matrix_path, genotype_path, locpath)):
            audit_failure(t, q, pairs, "missing_input_file")
            print(f"Genotype comparison inputs missing: {t} {q}: {matrix_path}, {genotype_path}, {locpath}", flush=True)
            continue
        wanted = {snp for snp, trait in pairs}
        locations = {}
        with locpath.open() as handle:
            for r in csv.DictReader(handle, delimiter="\t"):
                if r["snpid"] in wanted:
                    k = (chr_key(r["chr"]), int(r["pos"]))
                    if r["snpid"] in locations and locations[r["snpid"]] != k:
                        raise ValueError(f"Ambiguous SNP location: {r['snpid']}")
                    locations[r["snpid"]] = k
        jobs.append((t, q, pairs, matrix_path, scale, genotype_path, locations))

# Read the PLINK header once; reject ambiguous coordinates
alleles = defaultdict(set)
raw = ROOT/"QTL/plink/joint_all_chrs_matrixEQTL.raw"
targets = {k for job in jobs for k in job[-1].values()}
if jobs and raw.is_file():
    with raw.open() as handle:
        for token in handle.readline().split()[6:]:
            m = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([^:]+):([^_]+)_(.+)", token)
            if m:
                c, pos, ref, alt, counted = m.groups()
                k = (chr_key(c), int(pos))
                if k in targets:
                    alleles[k].add((ref, alt, counted))

for t, q, pairs, matrix_path, scale, genotype_path, locations in jobs:
    if not raw.is_file():
        audit_failure(t, q, pairs, "missing_plink_allele_header")
        continue
    subjects, phenotypes = selected_matrix(matrix_path, {trait for snp, trait in pairs})
    gs, genotypes = selected_matrix(genotype_path, {snp for snp, trait in pairs})
    if set(gs) != set(subjects):
        raise ValueError(f"Phenotype/genotype subject sets differ: {t} {q}")
    gi = {v: i for i, v in enumerate(gs)}
    order = [gi[v] for v in subjects]
    for snp, trait in pairs:
        r = dict(SNP=snp, tissue=t, qtl_type=q, trait=trait,
                 matrix_path=str(matrix_path), input_scale=scale,
                 threshold=THRESHOLDS[q], passes=False)
        aa = alleles.get(locations.get(snp), set())
        if trait not in phenotypes or snp not in genotypes:
            r["status"] = "missing_matrix_row"
        elif len(aa) != 1:
            r["status"] = "missing_or_ambiguous_alleles"
        else:
            ref, alt, counted = next(iter(aa))
            if counted not in (ref, alt) or "," in alt or ref == alt:
                r["status"] = "unsupported_alleles"
            else:
                x = phenotypes[trait]
                g = genotypes[snp][order]
                if np.any(np.isfinite(g) & ~np.isin(g, [0, 1, 2])):
                    raise ValueError(f"Non-hardcall genotype: {t} {q} {snp}")
                if scale != "zscore":
                    limit = 1 if scale == "fraction" else 100
                    if np.any(np.isfinite(x) & ((x < 0) | (x > limit))):
                        raise ValueError(f"Phenotype outside declared {scale} range: {t} {q} {trait}")
                altg = 2-g if counted == ref else g
                x = x * 100 if scale == "fraction" else x
                r.update(ref=ref, alt=alt, plink_counted_allele=counted,
                         units="zscore" if q == "eQTL" else "percentage_points")
                for value, label in [(0, "ref_ref"), (1, "het"), (2, "hom_alt")]:
                    v = x[(altg == value) & np.isfinite(x)]
                    r["n_"+label] = len(v)
                    r["mean_"+label] = float(v.mean()) if len(v) else np.nan
                r["delta_alt_minus_ref"] = r["mean_hom_alt"] - r["mean_ref_ref"]
                r["abs_delta"] = abs(r["delta_alt_minus_ref"])
                if min(r["n_ref_ref"], r["n_hom_alt"]) < MIN_HOMOZYGOTES:
                    r["status"] = "fewer_than_5_in_one_or_both_homozygote_groups"
                else:
                    r["status"] = "evaluated"
                    r["passes"] = bool(r["abs_delta"] >= THRESHOLDS[q])
        audit.append(r)

columns = ["SNP", "tissue", "qtl_type", "trait", "status", "passes", "matrix_path", "input_scale",
           "threshold", "ref", "alt", "plink_counted_allele", "units", "n_ref_ref", "n_het",
           "n_hom_alt", "mean_ref_ref", "mean_het", "mean_hom_alt", "delta_alt_minus_ref", "abs_delta"]
differences = pd.DataFrame(audit).reindex(columns=columns)
differences.to_csv(OUTDIR/"ALS_genotype_difference_audit.tsv", sep="\t", index=False, na_rep="NA")
evaluated_snps = set(differences.loc[differences.status.eq("evaluated"), "SNP"])
passing = differences[differences.passes.eq(True)]
passing_snps = set(passing.SNP)
out["large_genotype_difference"] = [snp in passing_snps for snp in out.index]
out["large_genotype_difference_status"] = [
    "PASS" if snp in passing_snps else "NO_PASS_AMONG_EVALUATED" if snp in evaluated_snps else "NOT_EVALUATED"
    for snp in out.index
]
pass_details = defaultdict(set)
for r in passing.to_dict("records"):
    pass_details[r["SNP"]].add(f"{r['qtl_type']}|{r['tissue']}|{r['trait']}")
out["large_genotype_difference_matches"] = [";".join(sorted(pass_details[snp])) for snp in out.index]
for q in QTLS:
    sub = differences[(differences.qtl_type == q) & differences.status.eq("evaluated")]
    best = sub.groupby("SNP").abs_delta.max()
    out[f"{q}_max_abs_homozygote_difference"] = out["SNP"].map(best)
print(f"Genotype difference audit: {OUTDIR/'ALS_genotype_difference_audit.tsv'}", flush=True)
print("Genotype comparison statuses:", differences.status.value_counts().to_dict(), flush=True)

# SCREEN annotation: all significant sQTL overlaps, without top-10 filtering.
screen_file = ROOT/"QTL/SCREEN_sQTL_PSI_ranking/SCREEN_overlaps.tsv"
if not screen_file.is_file():
    raise FileNotFoundError(f"Required SCREEN overlap table missing: {screen_file}. Run the SCREEN script first.")
screen = pd.read_csv(screen_file, sep="\t", dtype=str, keep_default_na=False)
required_screen = {"tissue", "junction_id", "snpid", "regulatory_class", "cCRE_id", "interval_bed"}
if not required_screen.issubset(screen.columns):
    raise ValueError(f"SCREEN table requires columns: {sorted(required_screen)}")
screen = screen[screen["tissue"].isin(TISSUES) & screen["regulatory_class"].isin(["PLS", "pELS", "dELS"])]
screen_sets = {k: defaultdict(set) for k in ["tissues", "classes", "cCRE_ids", "intervals_bed", "traits", "matches"]}
for r in screen.to_dict("records"):
    snp = clean_snp(r["snpid"])
    if snp not in out.index:
        continue
    for key, value in [("tissues", r["tissue"]), ("classes", r["regulatory_class"]),
                       ("cCRE_ids", r["cCRE_id"]), ("intervals_bed", r["interval_bed"]),
                       ("traits", r["junction_id"]),
                       ("matches", "|".join(r[k] for k in ["tissue", "junction_id", "regulatory_class", "cCRE_id", "interval_bed"]))]:
        if value:
            screen_sets[key][snp].add(value)
out["SCREEN_regulatory_overlap"] = [bool(screen_sets["matches"].get(snp)) for snp in out.index]
for key, values in screen_sets.items():
    out[f"SCREEN_{key}"] = [";".join(sorted(values.get(snp, set()))) for snp in out.index]

# ------------------------------------------------------------------
# 6. EVIDENCE COUNTS
#
# GWAS p<0.05 is now one evidence point.
# ------------------------------------------------------------------

score_cols=[
    "GWAS_sig",
    "ALS_relevant_gene",
    "NEUROLOGY_relevant_gene",
    "RNA_binding_protein",
    "large_genotype_difference",
    "SCREEN_regulatory_overlap",

    "eQTL_significant",
    "sQTL_significant",
    "cQTL_significant",

    "eQTL_SMR_significant",
    "sQTL_SMR_significant",
    "cQTL_SMR_significant",

    "eQTL_HEIDI_pass",
    "sQTL_HEIDI_pass",
    "cQTL_HEIDI_pass",

    "eQTL_coloc_significant",
    "sQTL_coloc_significant",
    "cQTL_coloc_significant",

    "eQTL_SuSiE_credible_set",
    "sQTL_SuSiE_credible_set",
    "cQTL_SuSiE_credible_set"
]

pd.DataFrame({"scored_category": score_cols, "points_if_true": 1}).to_csv(
    OUTDIR/"ALS_evidence_scoring_categories.tsv", sep="\t", index=False)
out["evidence_count"]=(
    out[
        score_cols
    ]
    .astype(bool)
    .astype(int)
    .sum(axis=1)
)

out["any_QTL"]=(
    out["eQTL_significant"]
    |
    out["sQTL_significant"]
    |
    out["cQTL_significant"]
)

out["any_SMR"]=(
    out["eQTL_SMR_significant"]
    |
    out["sQTL_SMR_significant"]
    |
    out["cQTL_SMR_significant"]
)

out["any_HEIDI_pass"]=(
    out["eQTL_HEIDI_pass"]
    |
    out["sQTL_HEIDI_pass"]
    |
    out["cQTL_HEIDI_pass"]
)

out["any_coloc"]=(
    out["eQTL_coloc_significant"]
    |
    out["sQTL_coloc_significant"]
    |
    out["cQTL_coloc_significant"]
)

out["any_SuSiE_credible_set"]=(
    out["eQTL_SuSiE_credible_set"]
    |
    out["sQTL_SuSiE_credible_set"]
    |
    out["cQTL_SuSiE_credible_set"]
)

out["any_evidence"]=(
    out["evidence_count"]>0
)

out["evidence_categories"] = [
    ";".join(c for c in score_cols if bool(row[c]))
    for row in out[score_cols].to_dict("records")
]

# ------------------------------------------------------------------
# 7. VALIDATION PRIORITY
#
# Transparent heuristic only.
#
# HIGH:
#   GWAS p<0.05 + QTL + coloc + SMR
#
# MEDIUM:
#   GWAS p<0.05 + QTL + (SMR or coloc)
#
# LOW:
#   GWAS p<0.05 + QTL
#
# GWAS_ONLY:
#   GWAS p<0.05 but no significant molecular QTL evidence
#
# MOLECULAR_ONLY:
#   molecular evidence exists but GWAS p>=0.05
# ------------------------------------------------------------------

out["validation_priority"]=np.select(
    [
        (
            out["GWAS_sig"]
            &
            out["any_QTL"]
            &
            out["any_SMR"]
            &
            out["any_coloc"]
        ),

        (
            out["GWAS_sig"]
            &
            out["any_QTL"]
            &
            (
                out["any_SMR"]
                |
                out["any_coloc"]
            )
        ),

        (
            out["GWAS_sig"]
            &
            out["any_QTL"]
        ),

        (
            out["GWAS_sig"]
        ),

        (
            ~out["GWAS_sig"]
            &
            (
                out["any_QTL"]
                |
                out["any_SMR"]
                |
                out["any_coloc"]
                |
                out["any_SuSiE_credible_set"]
            )
        )
    ],
    [
        "HIGH",
        "MEDIUM",
        "LOW",
        "GWAS_ONLY",
        "MOLECULAR_ONLY"
    ],
    default="NONE"
)

# ------------------------------------------------------------------
# 8. REORDER IMPORTANT COLUMNS
# ------------------------------------------------------------------

front_cols=[
    "SNP",

    "GWAS_chr",
    "GWAS_pos",
    "GWAS_effect_allele",
    "GWAS_other_allele",
    "GWAS_beta",
    "GWAS_se",
    "GWAS_p",
    "GWAS_sig",

    "candidate_gene_symbols",
    "candidate_gene_ids",
    "candidate_genes"
]

remaining=[
    c
    for c in out.columns
    if c not in front_cols
]

out=out[
    front_cols
    +
    remaining
]

# ------------------------------------------------------------------
# 9. OUTPUT
# ------------------------------------------------------------------

out=out.reset_index(
    drop=True
)

allfile=(
    OUTDIR
    /
    "ALS_GWAS_variant_validation_evidence.ALL.tsv.gz"
)

evfile=(
    OUTDIR
    /
    "ALS_GWAS_variant_validation_evidence.WITH_EVIDENCE.tsv"
)

rankfile=(
    OUTDIR
    /
    "ALS_GWAS_variant_validation_evidence.RANKED.tsv"
)

out.to_csv(
    allfile,
    sep="\t",
    index=False,
    compression="gzip"
)

evidence=out[
    out["any_evidence"]
].copy()

evidence.to_csv(
    evfile,
    sep="\t",
    index=False
)

priority_order={
    "HIGH":0,
    "MEDIUM":1,
    "LOW":2,
    "GWAS_ONLY":3,
    "MOLECULAR_ONLY":4,
    "NONE":5
}

ranked=evidence.copy()

ranked["_priority_order"]=(
    ranked[
        "validation_priority"
    ]
    .map(priority_order)
    .fillna(99)
)

ranked=ranked.sort_values(
    [
        "_priority_order",
        "evidence_count",
        "GWAS_p"
    ],
    ascending=[
        True,
        False,
        True
    ],
    na_position="last"
)

ranked=ranked.drop(
    columns=[
        "_priority_order"
    ]
)

ranked.insert(0, "priority_rank", np.arange(1, len(ranked)+1))
ranked.to_csv(rankfile, sep="\t", index=False)
# ------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------

print()
print("="*72)
print("DONE")
print("="*72)

print(
    f"GWAS SNPs                  : {len(out):,}"
)

print(
    f"GWAS p<0.05                : {int(out['GWAS_sig'].sum()):,}"
)

print(
    f"SNPs with QTL evidence     : {int(out['any_QTL'].sum()):,}"
)

print(
    f"SNPs with SMR evidence     : {int(out['any_SMR'].sum()):,}"
)

print(
    f"SNPs with HEIDI pass       : {int(out['any_HEIDI_pass'].sum()):,}"
)

print(
    f"SNPs with coloc            : {int(out['any_coloc'].sum()):,}"
)

print(
    f"SNPs in SuSiE credible set : {int(out['any_SuSiE_credible_set'].sum()):,}"
)

print(f"SNPs with SCREEN overlap    : {int(out['SCREEN_regulatory_overlap'].sum()):,}")
print(
    f"SNPs with candidate genes  : {(out['candidate_gene_ids']!='').sum():,}"
)

print(
    f"SNPs with any evidence     : {int(out['any_evidence'].sum()):,}"
)

print()
print("Validation priority:")

for level in [
    "HIGH",
    "MEDIUM",
    "LOW",
    "GWAS_ONLY",
    "MOLECULAR_ONLY",
    "NONE"
]:
    n=(
        out["validation_priority"]
        ==level
    ).sum()

    print(
        f"  {level:<15}: {n:,}"
    )

print()
print(
    f"ALL      : {allfile}"
)

print(
    f"EVIDENCE : {evfile}"
)

print(
    f"RANKED   : {rankfile}"
)

print("="*72)

PY
