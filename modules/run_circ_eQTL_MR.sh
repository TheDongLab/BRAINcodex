#!/usr/bin/env bash
#SBATCH --job-name=circ_eQTL_MR
#SBATCH --output=/home/zw529/donglab/data/target_ALS/MR/circ_eQTL_MR.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/MR/circ_eQTL_MR.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail
module --force purge
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
python3 -c 'import pandas,numpy' 2>/dev/null || { echo "ERROR: pandas/numpy unavailable"; exit 1; }

# Set EXPERIMENT=target_ALS or answer_ALS. TargetALS is the default.
EXPERIMENT="${EXPERIMENT:-target_ALS}"
[[ "$EXPERIMENT" == "target_ALS" || "$EXPERIMENT" == "answer_ALS" ]] || { echo "ERROR: EXPERIMENT must be target_ALS or answer_ALS" >&2; exit 1; }

BASE="$HOME/donglab/data/$EXPERIMENT"
GLOBAL="$BASE/MR/circ_eQTL_MR"
# Requires circRNA QTL effects already harmonized to GWAS effect alleles.
# Same input filename as the gene workflow; geneid or circ_id must hold circ IDs.
ROOT="${CIRC_MR_ROOT:-$HOME/donglab/data/GCST90027163/GWAS/circ_eQTL_harmonization}"

if [[ "$EXPERIMENT" == "target_ALS" ]]; then
    TISSUES=(Spinal_Cord Cortex Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)
else
    TISSUES=(iPSC_Motor_Neuron)
fi
[[ -z "${QTL_TISSUES:-}" ]] || read -r -a TISSUES <<< "$QTL_TISSUES"
[[ -d "$ROOT" ]] || { echo "ERROR: MR input root does not exist: $ROOT" >&2; exit 1; }
MIN_F=10
export BASE ROOT GLOBAL MIN_F TISSUE_LIST="${TISSUES[*]}"

python3 <<'PY'
import os, math
import numpy as np, pandas as pd

BASE,ROOT,GLOBAL=os.environ["BASE"],os.environ["ROOT"],os.environ["GLOBAL"]
TISSUES=os.environ["TISSUE_LIST"].split(); MIN_F=float(os.environ["MIN_F"])
os.makedirs(GLOBAL,exist_ok=True)

def bh(p):
    p=pd.to_numeric(p,errors="coerce"); out=pd.Series(np.nan,index=p.index,dtype=float); ok=p.notna()
    if not ok.any(): return out
    x=p[ok].clip(0,1); order=x.sort_values().index; vals=x.loc[order].values; m=len(vals); adj=np.empty(m); run=1.
    for i in range(m-1,-1,-1):
        run=min(run,vals[i]*m/(i+1)); adj[i]=min(run,1.)
    out.loc[order]=adj
    return out

def exp_safe(x):
    if pd.isna(x): return np.nan
    if x>700: return np.inf
    if x<-700: return 0.
    return math.exp(x)

def unique_join(x):
    return ";".join(sorted(set(map(str,x.dropna()))))

relationship_by_tissue={}; snp_by_tissue={}; summaries=[]

# ==============================================================
# 1. BUILD SNP×CIRCRNA EFFECT TABLE + UNIQUE-SNP STATISTICAL TABLE
# ==============================================================
for tissue in TISSUES:
    infile=f"{ROOT}/{tissue}/{tissue}_SNP_gene_MR_ready.tsv.gz"
    out=f"{BASE}/{tissue}/MR/circ_eQTL_MR"; os.makedirs(out,exist_ok=True)
    if not os.path.isfile(infile) or os.path.getsize(infile)==0: raise FileNotFoundError(f"Missing/empty: {infile}")

    print(f"\n{'='*70}\n{tissue}\n{'='*70}")
    d=pd.read_csv(infile,sep="\t",compression="gzip")
    # Preserve each BSJ as an exposure; never aggregate by host gene.
    if "circ_id" not in d.columns and "geneid" in d.columns:
        d=d.rename(columns={"geneid":"circ_id"})
    if "circ_id" in d.columns:
        valid=d["circ_id"].astype("string").str.fullmatch(r"[^:]+:[0-9]+-[0-9]+:[+-]",na=False)
        if not valid.all():
            raise ValueError(f"{tissue}: expected circRNA coordinate IDs, not host-gene/Ensembl IDs")
        if "geneid" in d.columns and not d["geneid"].astype(str).equals(d["circ_id"].astype(str)):
            raise ValueError(f"{tissue}: conflicting geneid and circ_id columns")


    req=["snpid","circ_id","beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]
    miss=set(req)-set(d.columns)
    if miss: raise ValueError(f"{tissue} input missing: {sorted(miss)}")

    for c in ["beta","se_qtl","p-value","gwas_beta_harmonized","gwas_se","gwas_p"]:
        d[c]=pd.to_numeric(d[c],errors="coerce")

    d=d.rename(columns={"beta":"beta_eqtl","p-value":"p_eqtl"})
    d["tissue"]=tissue
    d["F_statistic"]=(d["beta_eqtl"]/d["se_qtl"])**2

    good=(d["F_statistic"]>=MIN_F)&d["beta_eqtl"].notna()&(d["beta_eqtl"]!=0)&(d["se_qtl"]>0)&d["gwas_beta_harmonized"].notna()&(d["gwas_se"]>0)&d["gwas_p"].notna()
    d=d[good].copy().drop_duplicates(["snpid","circ_id"])

    # Wald ratio = effect-size annotation only
    d["beta_Wald_ratio"]=d["gwas_beta_harmonized"]/d["beta_eqtl"]
    d["se_Wald_first_order"]=d["gwas_se"]/d["beta_eqtl"].abs()
    d["se_Wald_delta"]=np.sqrt(
        d["gwas_se"]**2/d["beta_eqtl"]**2 +
        d["gwas_beta_harmonized"]**2*d["se_qtl"]**2/d["beta_eqtl"]**4
    )
    d["OR_Wald_ratio"]=d["beta_Wald_ratio"].apply(exp_safe)
    d["CI95_lower_beta"]=d["beta_Wald_ratio"]-1.96*d["se_Wald_first_order"]
    d["CI95_upper_beta"]=d["beta_Wald_ratio"]+1.96*d["se_Wald_first_order"]
    d["OR_CI95_lower"]=d["CI95_lower_beta"].apply(exp_safe)
    d["OR_CI95_upper"]=d["CI95_upper_beta"].apply(exp_safe)

    # UNIQUE SNP = statistical unit
    check=d.groupby("snpid").agg(
        n_beta=("gwas_beta_harmonized","nunique"),
        n_se=("gwas_se","nunique"),
        n_p=("gwas_p","nunique")
    )
    bad=check[(check.n_beta>1)|(check.n_se>1)|(check.n_p>1)]
    if len(bad): raise ValueError(f"{tissue}: {len(bad)} SNPs have inconsistent GWAS statistics across circRNA annotations")

    snp=d.sort_values(["snpid","p_eqtl"]).drop_duplicates("snpid").copy()
    keep=["tissue","snpid","gwas_beta_harmonized","gwas_se","gwas_p"]
    for c in ["gwas_chr","gwas_pos","harmonized_effect_allele","harmonized_other_allele"]:
        if c in snp.columns: keep.append(c)
    snp=snp[keep].copy()

    annot=d.groupby("snpid").agg(
        n_regulatory_targets=("circ_id","nunique"),
        regulatory_circRNAs=("circ_id",unique_join)
    ).reset_index()
    snp=snp.merge(annot,on="snpid",how="left",validate="one_to_one")

    snp["FDR_ALS_tissue"]=bh(snp["gwas_p"])
    snp=snp.sort_values(["FDR_ALS_tissue","gwas_p","snpid"])

    relationship_by_tissue[tissue]=d
    snp_by_tissue[tissue]=snp

    summaries.append({
        "tissue":tissue,
        "SNP_circRNA_relationships":len(d),
        "unique_SNPs":len(snp),
        "target_circRNAs":d["circ_id"].nunique(),
        "multi_target_SNPs":int((snp["n_regulatory_targets"]>=2).sum()),
        "ALS_nominal_P0.05_SNPs":int((snp["gwas_p"]<.05).sum()),
        "ALS_tissue_FDR0.05_SNPs":int((snp["FDR_ALS_tissue"]<.05).sum())
    })

    print(f"SNP→circRNA relationships : {len(d):,}")
    print(f"Unique SNPs            : {len(snp):,}")
    print(f"Target circRNAs           : {d['circ_id'].nunique():,}")
    print(f"Multi-target SNPs      : {(snp['n_regulatory_targets']>=2).sum():,}")
    print(f"ALS P<0.05 SNPs        : {(snp['gwas_p']<.05).sum():,}")
    print(f"Tissue FDR<0.05 SNPs  : {(snp['FDR_ALS_tissue']<.05).sum():,}")

# ==============================================================
# 2. GLOBAL FDR — EACH SNP COUNTED EXACTLY ONCE ACROSS ALL TISSUES
# ==============================================================
combined=pd.concat(snp_by_tissue.values(),ignore_index=True)

check=combined.groupby("snpid").agg(
    n_beta=("gwas_beta_harmonized","nunique"),
    n_se=("gwas_se","nunique"),
    n_p=("gwas_p","nunique")
)
bad=check[(check.n_beta>1)|(check.n_se>1)|(check.n_p>1)]
if len(bad): raise ValueError(f"Across tissues: {len(bad)} SNPs have inconsistent GWAS statistics")

global_snp=combined.sort_values("snpid").drop_duplicates("snpid")[["snpid","gwas_beta_harmonized","gwas_se","gwas_p"]].copy()
ga=combined.groupby("snpid").agg(
    n_tissues=("tissue","nunique"),
    tissues=("tissue",unique_join)
).reset_index()

all_rel=pd.concat(relationship_by_tissue.values(),ignore_index=True)
gg=all_rel.groupby("snpid").agg(
    n_unique_regulatory_circRNAs=("circ_id","nunique"),
    regulatory_circRNAs_all_tissues=("circ_id",unique_join)
).reset_index()

global_snp=global_snp.merge(ga,on="snpid",how="left").merge(gg,on="snpid",how="left")
global_snp["FDR_ALS_global"]=bh(global_snp["gwas_p"])
global_snp=global_snp.sort_values(["FDR_ALS_global","gwas_p","snpid"])
global_fdr=global_snp.set_index("snpid")["FDR_ALS_global"]

# ==============================================================
# 3. WRITE CONSOLIDATED circ-eQTL OUTPUTS
# ==============================================================
for tissue in TISSUES:
    out=f"{BASE}/{tissue}/MR/circ_eQTL_MR"

    snp=snp_by_tissue[tissue].copy()
    snp["FDR_ALS_global"]=snp["snpid"].map(global_fdr)
    snp=snp.sort_values(["FDR_ALS_global","FDR_ALS_tissue","gwas_p","snpid"])

    rel=relationship_by_tissue[tissue].copy()
    rel=rel.merge(
        snp[["snpid","n_regulatory_targets","FDR_ALS_tissue","FDR_ALS_global"]],
        on="snpid",how="left",validate="many_to_one"
    )
    rel=rel.sort_values(["FDR_ALS_global","gwas_p","snpid","circ_id"])

    # 1 row per SNP
    snp.to_csv(f"{out}/{tissue}_circ_eQTL_SNP_results.tsv.gz",sep="\t",index=False,compression="gzip")

    # 1 row per SNP×circRNA
    rel.to_csv(f"{out}/{tissue}_circ_eQTL_SNP_circ_effects.tsv.gz",sep="\t",index=False,compression="gzip")

global_snp.to_csv(f"{GLOBAL}/all_tissues_circ_eQTL_SNP_results.tsv.gz",sep="\t",index=False,compression="gzip")

summary=pd.DataFrame(summaries)
summary.to_csv(f"{GLOBAL}/circ_eQTL_MR_summary_all_tissues.tsv",sep="\t",index=False)

print(f"\n{'='*70}\nALL-TISSUE circ-eQTL SNP-FIRST SUMMARY\n{'='*70}")
print(summary.to_string(index=False))
print(f"\nUnique SNPs globally          : {len(global_snp):,}")
print(f"ALS P<0.05 globally           : {(global_snp['gwas_p']<.05).sum():,}")
print(f"Global SNP-level FDR<0.05     : {(global_snp['FDR_ALS_global']<.05).sum():,}")
PY

echo
echo "circ-eQTL SNP-first analysis complete."
column -t -s $'\t' "$GLOBAL/circ_eQTL_MR_summary_all_tissues.tsv"
