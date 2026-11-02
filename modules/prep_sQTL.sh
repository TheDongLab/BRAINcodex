#!/bin/bash
#SBATCH --job-name=prep_sQTL
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/%x_%j.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=649G

set -euo pipefail
module load BCFtools

[[ $# -ge 1 ]] || { echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch prep_sQTL.sh TISSUE" >&2; exit 1; }
EXPERIMENT="${EXPERIMENT:-target_ALS}"
TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')
[[ "$EXPERIMENT" == target_ALS || "$EXPERIMENT" == answer_ALS ]] || { echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2; exit 1; }

BASE="/home/zw529/donglab/data/$EXPERIMENT"
QTL_DIR="$BASE/QTL"
PLINK="$QTL_DIR/plink"
SPLICING_RAW="$QTL_DIR/splicing_matrix.txt"
RAW="$PLINK/joint_all_chrs_matrixEQTL.raw"
BIM="$PLINK/joint_all_chrs_filtered_bed.bim"
PCA="$PLINK/joint_pca.eigenvec"
MAP_FILE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/ncbi_to_ucsc.txt"
DBSNP_VCF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/GCF_000001405.40.gz"

if [[ "$EXPERIMENT" == target_ALS ]]; then
    OUTDIR="$BASE/$TISSUE_DIR/sQTL"
    METADATA_CSV="$BASE/targetALS_rnaseq_metadata.csv"
    BREAKDOWN_TSV="$QTL_DIR/patient_tissue_breakdown.tsv"
    COV_FILE="$QTL_DIR/covariates.tsv"
    SEX_MISMATCH_FILE="$BASE/$TISSUE_DIR/eQTL/potential_sex_mismatches.txt"
    SPLICING_LOC_SRC="$OUTDIR/tmp_splicing_events_hg38.bed"
else
    [[ "$TISSUE_DIR" == iPSC_Motor_Neuron ]] || { echo "[ERROR] AnswerALS branch supports iPSC_Motor_Neuron only." >&2; exit 1; }
    OUTDIR="$BASE/iPSC_Motor_Neuron/sQTL"
    METADATA_CSV="$BASE/aals_dataportal_datatable.csv"
    BREAKDOWN_TSV=""
    COV_FILE=""
    SEX_MISMATCH_FILE="$BASE/iPSC_Motor_Neuron/eQTL/potential_sex_mismatches.txt"
    SPLICING_LOC_SRC="$OUTDIR/tmp_splicing_events_hg38.bed"
fi
mkdir -p "$OUTDIR"

for f in "$SPLICING_RAW" "$RAW" "$BIM" "$PCA" "$MAP_FILE" "$DBSNP_VCF" "$METADATA_CSV"; do
    [[ -s "$f" ]] || { echo "[ERROR] Missing or empty input: $f" >&2; exit 1; }
done

echo "============================================"
echo "  $EXPERIMENT sQTL Prep for $TISSUE"
echo "============================================"

TMP_TARGETS="$OUTDIR/tmp_targets_ncbi.txt"
TMP_RSID_MAP="$OUTDIR/tmp_coord_to_rsid.txt"

awk -v map_file="$MAP_FILE" '
BEGIN {
    OFS="\t"
    while ((getline < map_file) > 0) map[toupper($2)] = toupper($1)
}
{
    ucsc = ($1 ~ /^[Cc][Hh][Rr]/) ? toupper($1) : "CHR"toupper($1)
    ncbi = (ucsc in map) ? map[ucsc] : ucsc
    print ncbi, $4, $4
}' "$BIM" > "$TMP_TARGETS"

bcftools query -T "$TMP_TARGETS" -f '%CHROM:%POS\t%ID\n' "$DBSNP_VCF" > "$TMP_RSID_MAP"
rm -f "$TMP_TARGETS"

echo "[0] Generating splicing-event location source..."
python3 << EOF
import pandas as pd
matrix = pd.read_csv("$SPLICING_RAW", sep="\\t", usecols=[0])
with open("$SPLICING_LOC_SRC", "w") as out:
    for event in matrix.iloc[:, 0].astype(str):
        try:
            chrom, strand, interval = event.split(":")
            start, end = interval.split("-")
            out.write(f"{chrom}\\t{start}\\t{end}\\t{event}\\n")
        except ValueError:
            continue
EOF

echo "[1] Intersecting metadata, genotypes, and splicing matrix..."
python3 << EOF
import pandas as pd, numpy as np, re, os, sys
from collections import Counter

experiment, tissue = "$EXPERIMENT", "$TISSUE"
raw = pd.read_csv("$RAW", sep=r"\\s+", usecols=lambda c: c not in ["PAT","MAT","SEX","PHENOTYPE"])
pca = pd.read_csv("$PCA", sep=r"\\s+").rename(columns={"#IID":"IID"})
headers = pd.read_csv("$SPLICING_RAW", sep="\\t", nrows=0).columns.tolist()
mismatch_file = "$SEX_MISMATCH_FILE"
mismatches = set()
if os.path.isfile(mismatch_file):
    mismatches = {x.strip() for x in open(mismatch_file) if x.strip()}

def read_rsid_map():
    d = {}
    for line in open("$TMP_RSID_MAP"):
        fields = line.rstrip("\\n").split("\\t")
        if len(fields) == 2 and fields[1] not in ("", "."):
            d[fields[0].upper()] = fields[1]
    return d

coord_to_rsid = read_rsid_map()
ucsc_to_ncbi = {}
for line in open("$MAP_FILE"):
    if line.strip():
        ncbi, ucsc = line.strip().split("\\t")
        ucsc_to_ncbi[ucsc.upper()] = ncbi.upper()

if experiment == "target_ALS":
    breakdown = pd.read_csv("$BREAKDOWN_TSV", sep="\\t") if tissue not in ("Cortex","Spinal_Cord") else None
    meta = pd.read_csv("$METADATA_CSV", low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower()
    meta = meta[~meta["externalsampleid"].isin(mismatches)].copy()
    patterns = {
        "Cortex": r"(?:^|[\\W_])cortex(?:$|[\\W_])",
        "Spinal_Cord": r"(?:^|[\\W_])spinal[\\W_]*cord(?:$|[\\W_])",
        "Motor_Cortex": "Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4",
        "Cervical_Spinal_Cord": "Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical",
        "Lumbar_Spinal_Cord": "Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral",
        "Thoracic_Spinal_Cord": "Thoracic Spinal Cord|Thoracic",
        "Frontal_Cortex": "Frontal Cortex|Frontal",
        "Cerebellum": "Cerebellum"
    }
    pattern = patterns.get(tissue, tissue)
    if tissue in ("Cortex","Spinal_Cord"):
        mt = meta[meta["tissue"].astype("string").str.contains(pattern, case=False, na=False)].copy()
    else:
        approved = breakdown[breakdown["tissues"].str.contains(tissue, na=False)]["subject_id"].unique()
        mf = meta[meta["externalsubjectid"].isin(approved)].copy()
        mt = mf[mf.apply(lambda r: r.astype(str).str.contains(pattern, case=False).any(), axis=1)].copy()
    def rin(r):
        try: return max(float(r.iloc[16]) if pd.notnull(r.iloc[16]) else 0, float(r.iloc[17]) if pd.notnull(r.iloc[17]) else 0)
        except Exception: return 0
    mt["RIN_score"] = mt.apply(rin, axis=1)
    mt = mt[(mt["post_mortem_interval_in_hours"] <= 40) & (mt["RIN_score"] >= 3)]
    meta_unique = mt.sort_values("RIN_score", ascending=False).drop_duplicates("externalsubjectid").copy()
    meta_unique["donor_id"] = meta_unique["externalsubjectid"].astype(str)
    meta_unique["matrix_id"] = meta_unique["externalsampleid"].astype(str).str.replace("-", "_", regex=False)
else:
    if tissue != "iPSC_Motor_Neuron":
        raise SystemExit("AnswerALS branch supports iPSC_Motor_Neuron only.")
    meta = pd.read_csv("$METADATA_CSV", low_memory=False, encoding="utf-8-sig")
    id_col = meta.columns[0]
    sex_col = next((c for c in meta.columns if c.strip().lower() == "sex"), None)
    group_col = next((c for c in meta.columns if c.strip().lower() in ("subject_group","subject group")), None)
    if sex_col is None or group_col is None:
        raise SystemExit("AnswerALS metadata must contain Sex and SUBJECT_GROUP columns.")
    meta = meta.rename(columns={id_col:"donor_id", sex_col:"sex", group_col:"subject_group"})
    meta["donor_id"] = meta["donor_id"].astype(str).str.strip()
    meta = meta[(meta["donor_id"] != "") & (~meta["donor_id"].isin(mismatches))]
    meta_unique = meta.drop_duplicates("donor_id").copy()
    meta_unique["matrix_id"] = meta_unique["donor_id"]

# Map donors to matrix columns and verify their per-sample PSI files for AnswerALS.
common = sorted(set(meta_unique["donor_id"].astype(str)) & set(raw["IID"].astype(str)))
sub_to_matrix, final = {}, []
junction_counts = Counter()

if experiment == "target_ALS":
    encodings = ["utf-8","latin1","iso-8859-1","cp1252"]
    for donor in common:
        row = meta_unique[meta_unique["donor_id"].astype(str) == donor].iloc[0]
        prefix = row["matrix_id"]
        matches = [c for c in headers if c.startswith(prefix)]
        if matches:
            final.append(donor); sub_to_matrix[donor] = matches[0]
            raw_sample = row["externalsampleid"].replace("-", "_")
            paths = []
            for name in os.listdir("/home/zw529/donglab/data/target_ALS"):
                path = os.path.join("/home/zw529/donglab/data/target_ALS", name, "RNAseq", "Processed", raw_sample, "leafcutter", "psi", f"{raw_sample}.leafcutter.PSI.tsv")
                if os.path.exists(path): paths.append(path)
            if paths:
                try:
                    psi = pd.read_csv(paths[0], sep=r"\\s+", usecols=["chrom","strand","start","end","reads"])
                    ids = psi["chrom"].astype(str)+":"+psi["strand"].astype(str)+":"+psi["start"].astype(str)+"-"+psi["end"].astype(str)
                    junction_counts.update(ids[psi["reads"] >= 5].unique())
                except Exception: pass
    if not final:
        raise SystemExit("No TargetALS samples matched metadata, genotype, matrix, and PSI files.")
    min_pass = max(1, int(len(final)*0.10))
    valid_junctions = {j for j,c in junction_counts.items() if c >= min_pass}
else:
    for donor in common:
        if donor in headers:
            final.append(donor); sub_to_matrix[donor] = donor
    if not final:
        raise SystemExit("No AnswerALS donors matched metadata, genotype, and splicing matrix.")
    valid_junctions = None

target_cols = ["junction_id"] + [sub_to_matrix[s] for s in final]
if "junction_id" not in headers:
    raise SystemExit("Splicing matrix must have a junction_id column.")
splicing = pd.read_csv("$SPLICING_RAW", sep="\\t", usecols=target_cols, index_col=0, na_values=["NA"])
splicing = splicing[[sub_to_matrix[s] for s in final]]
splicing.columns = final
if valid_junctions is not None:
    splicing = splicing.loc[splicing.index.isin(valid_junctions)]
splicing = splicing.loc[splicing.notna().mean(axis=1) >= 0.80]
splicing = splicing.loc[splicing.var(axis=1, skipna=True) > 1e-10]

snps = raw[raw["IID"].astype(str).isin(final)].set_index("IID").drop(columns=["FID"]).loc[final].T
numeric = snps.apply(pd.to_numeric, errors="coerce")
called = numeric.notna().sum(axis=1)
af = numeric.sum(axis=1, skipna=True)/(2*called)
maf = np.minimum(af,1-af)
snps = snps.loc[(called > 0) & (maf >= 0.05)]
def groups_pass(row):
    counts = pd.to_numeric(row, errors="coerce").dropna().value_counts()
    return len(counts) >= 2 and counts.min() >= 5
snps = snps.loc[snps.apply(groups_pass, axis=1)]

def to_rsid(v):
    clean = v.rsplit("_",1)[0] if "_" in v else v
    parts = str(clean).split(":")
    if len(parts) < 2: return clean
    chrom,pos = parts[0],parts[1]
    ucsc = chrom if chrom.upper().startswith("CHR") else "chr"+chrom
    ncbi = ucsc_to_ncbi.get(ucsc.upper(),ucsc.upper())
    return coord_to_rsid.get(f"{ncbi}:{pos}".upper(),clean)
snps.index = [to_rsid(str(x)) for x in snps.index]
snps = snps.loc[~snps.index.duplicated(keep="first")]

if experiment == "target_ALS":
    cov_full = pd.read_csv("$COV_FILE", sep="\\t")
    cov = cov_full[["externalsampleid"]+[c for c in cov_full.columns if c not in meta.columns]].drop_duplicates("externalsampleid")
    aligned = meta_unique[meta_unique["donor_id"].isin(final)].copy()
    merged = pd.merge(pd.merge(aligned,cov,on="externalsampleid"),pca,left_on="donor_id",right_on="IID")
    merged = merged.sort_values("RIN_score",ascending=False).drop_duplicates("donor_id")
    merged["sex_bin"] = merged["sex"].astype(str).str.lower().map({"male":1,"female":0})
    merged["is_als"] = merged["subject_group"].apply(lambda x: 1 if "ALS" in str(x) else 0)
    covars = merged.set_index("donor_id").reindex(final)[["sex_bin","age_at_death","PC1","PC2","PC3","PC4","PC5","is_als"]].T
else:
    aligned = meta_unique[meta_unique["donor_id"].isin(final)].copy()
    merged = pd.merge(aligned,pca,left_on="donor_id",right_on="IID",how="inner")
    if len(merged) != len(final):
        raise SystemExit("Some AnswerALS donors lack genotype PCA rows.")
    merged["sex_bin"] = merged["sex"].astype(str).str.lower().map({"male":1,"female":0})
    merged["is_als"] = merged["subject_group"].astype(str).str.contains(r"\\bALS\\b",case=False,na=False).astype(int)
    covars = merged.set_index("donor_id").reindex(final)[["sex_bin","PC1","PC2","PC3","PC4","PC5","is_als"]].T

splicing.to_csv("$OUTDIR/splicing_${TISSUE_DIR}.txt",sep="\\t",index=True,index_label="geneid",na_rep="NA")
snps.to_csv("$OUTDIR/snp_${TISSUE_DIR}.txt",sep="\\t",index=True,index_label="snpid",na_rep="NA")
covars.to_csv("$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt",sep="\\t",index=True,index_label="id",na_rep="NA")
print(f"SUCCESS: Generated filtered datasets for {len(final)} samples; {len(splicing)} junctions, {len(snps)} SNPs.")
EOF

echo "Generating location files..."
printf "snpid\tchr\tpos\n" > "$OUTDIR/snp_location.txt"
awk -v ncbi_map="$MAP_FILE" -v rsid_map="$TMP_RSID_MAP" '
BEGIN {
    OFS="\t"
    while ((getline < ncbi_map) > 0) { n_map[toupper($2)] = toupper($1); inv_map[toupper($1)] = $2 }
    while ((getline < rsid_map) > 0) { r_map[toupper($1)] = $2 }
}
{
    split($2,p,":"); chrom=p[1]; pos=p[2]
    ucsc=(chrom ~ /^[Cc][Hh][Rr]/)?toupper(chrom):"CHR"toupper(chrom)
    ncbi=(ucsc in n_map)?n_map[ucsc]:ucsc
    id=(toupper(ncbi":"pos) in r_map)?r_map[toupper(ncbi":"pos)]:ncbi":"pos
    final_chr=(toupper(ncbi) in inv_map)?inv_map[toupper(ncbi)]:chrom
    if (!seen[id]++) print id,final_chr,pos
}' "$BIM" >> "$OUTDIR/snp_location.txt"

printf "geneid\tchr\tleft\tright\n" > "$OUTDIR/splicing_location.txt"
awk 'BEGIN{OFS="\t"} {if (!seen[$4]++) print $4,$1,$2,$3}' "$SPLICING_LOC_SRC" >> "$OUTDIR/splicing_location.txt"

rm -f "$SPLICING_LOC_SRC" "$TMP_RSID_MAP"
echo "Location files complete for $TISSUE."
