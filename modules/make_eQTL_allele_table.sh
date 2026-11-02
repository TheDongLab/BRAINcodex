#!/usr/bin/env bash
#SBATCH --job-name=make_eQTL_alleles
#SBATCH --output=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.out
#SBATCH --error=/home/zw529/donglab/data/GCST90027163/GWAS/make_eQTL_alleles.err
#SBATCH --time=2:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
command -v python3 >/dev/null || { echo "ERROR: python3 not found"; exit 1; }
python3 -c 'import pandas' 2>/dev/null || { echo "ERROR: pandas unavailable"; exit 1; }

# Set EXPERIMENT=target_ALS or answer_ALS; TargetALS remains the default.
EXPERIMENT="${EXPERIMENT:-target_ALS}"
[[ "$EXPERIMENT" == "target_ALS" || "$EXPERIMENT" == "answer_ALS" ]] || { echo "ERROR: EXPERIMENT must be target_ALS or answer_ALS" >&2; exit 1; }
BASE="$HOME/donglab/data/$EXPERIMENT"; PLINK_DIR="$BASE/QTL/plink"
BFILE="$PLINK_DIR/joint_all_chrs_filtered_bed"; RAW="$PLINK_DIR/joint_all_chrs_matrixEQTL.raw"; BIM="$BFILE.bim"
GWAS="$HOME/donglab/references/GWAS/ALS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
if [[ "$EXPERIMENT" == "target_ALS" ]]; then
    OUTROOT="$HOME/donglab/data/GCST90027163/GWAS/eQTL_harmonization"
    TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)
else
    OUTROOT="$BASE/MR/eQTL_harmonization"
    TISSUES=(iPSC_Motor_Neuron)
fi
SUMMARY="$OUTROOT/all_tissues_SNP_first_summary.tsv"
MIN_F=10

mkdir -p "$OUTROOT"
for f in "$RAW" "$BIM" "$GWAS"; do [[ -s "$f" ]] || { echo "ERROR: missing/empty: $f"; exit 1; }; done
export BASE RAW BIM GWAS OUTROOT SUMMARY MIN_F TISSUE_LIST="${TISSUES[*]}"

python3 <<'PY'
import os
import pandas as pd

BASE,RAW,BIM,GWAS,OUTROOT,SUMMARY=[os.environ[x] for x in ["BASE","RAW","BIM","GWAS","OUTROOT","SUMMARY"]]
TISSUES=os.environ["TISSUE_LIST"].split(); MIN_F=float(os.environ["MIN_F"])
META={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}; COMP={"A":"T","T":"A","C":"G","G":"C"}

def norm_chr(s): return s.astype(str).str.replace(r"^chr","",regex=True,case=False).str.upper()
def pct(n,d): return 100*n/d if d else 0
def req(df,cols,label):
    m=set(cols)-set(df.columns)
    if m: raise ValueError(f"{label} missing columns: {sorted(m)}")

def qtl_stats(d):
    d=d.copy()
    for c in ["beta","t-stat","p-value","FDR"]: d[c]=pd.to_numeric(d[c],errors="coerce")
    d["se_qtl"]=(d["beta"]/d["t-stat"]).abs()
    d["F_statistic"]=(d["beta"]/d["se_qtl"])**2
    return d

def harmonize(qtl,alleles,gwas):
    x=qtl.merge(alleles,on="snpid",how="left",validate="many_to_one")
    x=x.merge(gwas,on="snpid",how="left",validate="many_to_one")

    x["found_in_gwas"]=x["gwas_chr"].notna()
    x["allele_coordinate_match"]=(norm_chr(x["qtl_chr"])==norm_chr(x["allele_chr"]))&(x["qtl_pos"].astype(str)==x["allele_pos"].astype(str))
    x["gwas_coordinate_match"]=x["found_in_gwas"]&(norm_chr(x["qtl_chr"])==norm_chr(x["gwas_chr"]))&(x["qtl_pos"].astype(str)==x["gwas_pos"].astype(str))

    for c in ["effect_allele","other_allele","gwas_effect_allele","gwas_other_allele"]: x[c]=x[c].str.upper()
    x["palindromic"]=(x["effect_allele"].fillna("")+x["other_allele"].fillna("")).isin(["AT","TA","CG","GC"])

    def rel(r):
        if not r["found_in_gwas"]: return "not_found_in_GWAS"
        qe,qo,ge,go=r["effect_allele"],r["other_allele"],r["gwas_effect_allele"],r["gwas_other_allele"]
        if pd.isna(qe) or pd.isna(qo): return "missing_QTL_alleles"
        if pd.isna(ge) or pd.isna(go): return "missing_GWAS_alleles"
        if qe==ge and qo==go: return "exact"
        if qe==go and qo==ge: return "swap_beta"
        if qe==COMP.get(ge) and qo==COMP.get(go): return "strand"
        if qe==COMP.get(go) and qo==COMP.get(ge): return "strand_and_swap_beta"
        return "incompatible"

    x["harmonization_relation"]=x.apply(rel,axis=1)
    x["beta_flip_required"]=x["harmonization_relation"].isin(["swap_beta","strand_and_swap_beta"])
    x["strand_flip_required"]=x["harmonization_relation"].isin(["strand","strand_and_swap_beta"])
    x["gwas_beta_harmonized"]=x["gwas_beta_original"]
    x.loc[x["beta_flip_required"],"gwas_beta_harmonized"]=-x.loc[x["beta_flip_required"],"gwas_beta_original"]

    x["harmonized_effect_allele"]=x["effect_allele"]; x["harmonized_other_allele"]=x["other_allele"]; x["rejection_reason"]=""
    x.loc[x["harmonization_relation"]=="not_found_in_GWAS","rejection_reason"]="not_found_in_GWAS"
    x.loc[x["harmonization_relation"]=="missing_QTL_alleles","rejection_reason"]="missing_QTL_alleles"
    x.loc[x["harmonization_relation"]=="missing_GWAS_alleles","rejection_reason"]="missing_GWAS_alleles"
    x.loc[x["harmonization_relation"]=="incompatible","rejection_reason"]="incompatible_alleles"
    x.loc[(x["rejection_reason"]=="")&x["palindromic"],"rejection_reason"]="palindromic_without_frequency_check"
    x.loc[(x["rejection_reason"]=="")&~x["allele_coordinate_match"],"rejection_reason"]="QTL_allele_coordinate_mismatch"
    x.loc[(x["rejection_reason"]=="")&~x["gwas_coordinate_match"],"rejection_reason"]="GRCh38_coordinate_mismatch"
    x.loc[(x["rejection_reason"]=="")&(x["allele_status"]!="ok"),"rejection_reason"]="QTL_allele_status_not_ok"
    x.loc[(x["rejection_reason"]=="")&(x["beta"].isna()|x["se_qtl"].isna()|x["gwas_beta_harmonized"].isna()|x["gwas_se"].isna()),"rejection_reason"]="missing_beta_or_standard_error"
    x.loc[(x["rejection_reason"]=="")&((x["se_qtl"]<=0)|(x["gwas_se"]<=0)),"rejection_reason"]="nonpositive_standard_error"

    x["harmonization_keep"]=x["rejection_reason"]==""
    x["strong_instrument"]=x["F_statistic"]>=MIN_F
    x["MR_ready"]=x["harmonization_keep"]&x["strong_instrument"]
    return x

# ------------------------------------------------------------------
# Recover MatrixEQTL counted/effect allele from PLINK RAW header
# ------------------------------------------------------------------
print("Reading RAW header...")
with open(RAW) as f: cols=f.readline().split()

records=[]
for c in [x for x in cols if x not in META]:
    if "_" not in c: records.append((c,None,"missing_counted_allele_suffix"))
    else:
        vid,a=c.rsplit("_",1); records.append((vid,a.upper(),"ok"))

raw_meta=pd.DataFrame(records,columns=["raw_variant_id","effect_allele","raw_parse_status"])
print(f"RAW genotype variants: {len(raw_meta):,}")

# ------------------------------------------------------------------
# Recover opposite allele and GRCh38 coordinate from BIM
# ------------------------------------------------------------------
print("Reading BIM...")
bim=pd.read_csv(BIM,sep=r"\s+",header=None,names=["bim_chr","bim_id","cm","bim_pos","bim_a1","bim_a2"],dtype=str)
bim["bim_a1"]=bim["bim_a1"].str.upper(); bim["bim_a2"]=bim["bim_a2"].str.upper()

vm=raw_meta.merge(bim,left_on="raw_variant_id",right_on="bim_id",how="left",validate="one_to_one")
vm["other_allele"]=vm.apply(lambda r:r["bim_a2"] if r["effect_allele"]==r["bim_a1"] else r["bim_a1"] if r["effect_allele"]==r["bim_a2"] else None,axis=1)
vm["allele_status"]=vm.apply(lambda r:"ok" if r["effect_allele"] in {r["bim_a1"],r["bim_a2"]} else "missing_bim_match" if pd.isna(r["bim_id"]) else "counted_allele_not_in_bim",axis=1)
vm["key"]=norm_chr(vm["bim_chr"])+":"+vm["bim_pos"].astype(str)

# ------------------------------------------------------------------
# External ALS GWAS
# ------------------------------------------------------------------
print("Reading ALS GWAS...")
gc=["hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"]
gwas=pd.read_csv(GWAS,sep="\t",compression="gzip",usecols=gc,dtype={"hm_rsid":str,"hm_chrom":str,"hm_pos":str,"hm_other_allele":str,"hm_effect_allele":str})
gwas=gwas.rename(columns={"hm_rsid":"snpid","hm_chrom":"gwas_chr","hm_pos":"gwas_pos","hm_other_allele":"gwas_other_allele","hm_effect_allele":"gwas_effect_allele","hm_beta":"gwas_beta_original","standard_error":"gwas_se","p_value":"gwas_p"})
for c in ["gwas_effect_allele","gwas_other_allele"]: gwas[c]=gwas[c].str.upper()
for c in ["gwas_beta_original","gwas_se","gwas_p"]: gwas[c]=pd.to_numeric(gwas[c],errors="coerce")
gwas=gwas.drop_duplicates("snpid")

summary=[]

# ==================================================================
# TISSUE-SPECIFIC SNP-FIRST TABLES
# ==================================================================
for tissue in TISSUES:
    E=f"{BASE}/{tissue}/eQTL"; R=f"{E}/results"; O=f"{OUTROOT}/{tissue}"; os.makedirs(O,exist_ok=True)
    fdr_path=f"{R}/{tissue}_eQTL.FDR0.05.txt"; loc_path=f"{E}/snp_location.txt"
    for f in [fdr_path,loc_path]:
        if not os.path.isfile(f) or os.path.getsize(f)==0: raise FileNotFoundError(f"Missing/empty: {f}")

    print(f"\n{'='*70}\n{tissue}\n{'='*70}")

    loc=pd.read_csv(loc_path,sep="\t",dtype=str)
    req(loc,["snpid","chr","pos"],f"{tissue} locations")
    loc=loc.drop_duplicates("snpid"); loc["key"]=norm_chr(loc["chr"])+":"+loc["pos"].astype(str)

    qtl=pd.read_csv(fdr_path,sep="\t",dtype=str)
    req(qtl,["geneid","snpid","beta","t-stat","p-value","FDR"],f"{tissue} FDR0.05")
    qtl=qtl.merge(loc[["snpid","chr","pos"]],on="snpid",how="left",validate="many_to_one").rename(columns={"chr":"qtl_chr","pos":"qtl_pos"})
    qtl=qtl_stats(qtl).sort_values(["snpid","geneid","p-value"]).drop_duplicates(["snpid","geneid"])
    qtl.insert(0,"tissue",tissue)

    # Alleles are SNP-level, independent of which gene the SNP regulates
    needed=set(qtl["snpid"].dropna())
    alleles=loc[loc["snpid"].isin(needed)].merge(
        vm[["key","effect_allele","other_allele","bim_a1","bim_a2","raw_variant_id","allele_status"]],
        on="key",how="left"
    )
    alleles=alleles.rename(columns={"chr":"allele_chr","pos":"allele_pos","bim_a1":"bim_A1","bim_a2":"bim_A2"})
    alleles=alleles[["snpid","allele_chr","allele_pos","effect_allele","other_allele","bim_A1","bim_A2","raw_variant_id","allele_status"]].drop_duplicates("snpid")

    # Every row is an SNP x gene regulatory relationship
    x=harmonize(qtl,alleles,gwas)
    ready=x[x["MR_ready"]].copy()

    # Number of regulatory targets is annotation only; no gene aggregation occurs
    n_targets=ready.groupby("snpid")["geneid"].nunique()
    ready["n_regulatory_targets"]=ready["snpid"].map(n_targets).astype(int)

    # Pure SNP-level GWAS/genotype master: exactly one row per SNP
    snp_cols=["tissue","snpid","qtl_chr","qtl_pos","raw_variant_id",
              "effect_allele","other_allele","gwas_chr","gwas_pos",
              "gwas_effect_allele","gwas_other_allele","gwas_beta_original",
              "gwas_beta_harmonized","gwas_se","gwas_p",
              "harmonization_relation","beta_flip_required","strand_flip_required"]
    snp_master=ready[snp_cols].drop_duplicates("snpid").copy()
    snp_master["n_regulatory_targets"]=snp_master["snpid"].map(n_targets).astype(int)

    x.to_csv(f"{O}/{tissue}_SNP_gene_all.tsv.gz",sep="\t",index=False,compression="gzip")
    ready.to_csv(f"{O}/{tissue}_SNP_gene_MR_ready.tsv.gz",sep="\t",index=False,compression="gzip")
    snp_master.to_csv(f"{O}/{tissue}_SNP_GWAS_master.tsv.gz",sep="\t",index=False,compression="gzip")
    x[~x["MR_ready"]].to_csv(f"{O}/{tissue}_SNP_gene_rejected.tsv.gz",sep="\t",index=False,compression="gzip")

    summary.append({
        "tissue":tissue,
        "SNP_gene_FDR0.05_rows":len(qtl),
        "unique_eQTL_SNPs":qtl["snpid"].nunique(),
        "unique_target_genes":qtl["geneid"].nunique(),
        "harmonized_rows":int(x["harmonization_keep"].sum()),
        "strong_F10_rows":int((x["harmonization_keep"]&x["strong_instrument"]).sum()),
        "MR_ready_SNP_gene_rows":len(ready),
        "MR_ready_unique_SNPs":ready["snpid"].nunique(),
        "MR_ready_target_genes":ready["geneid"].nunique(),
        "multi_target_SNPs":int((n_targets>=2).sum()),
        "minimum_F":MIN_F
    })

    print(f"FDR-significant SNP×gene relationships : {len(qtl):,}")
    print(f"Unique eQTL SNPs                      : {qtl['snpid'].nunique():,}")
    print(f"MR-ready SNP×gene relationships       : {len(ready):,}")
    print(f"MR-ready unique SNPs                  : {ready['snpid'].nunique():,}")
    print(f"SNPs regulating >=2 genes             : {(n_targets>=2).sum():,}")
    print("No LD pruning or gene-level MR aggregation performed.")

pd.DataFrame(summary).to_csv(SUMMARY,sep="\t",index=False)
print("\nSNP-FIRST SUMMARY")
print(pd.DataFrame(summary).to_string(index=False))
PY

echo
echo "=== SNP-FIRST SUMMARY ==="
column -t -s $'\t' "$SUMMARY"
