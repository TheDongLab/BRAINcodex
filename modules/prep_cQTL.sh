#!/bin/bash
#SBATCH --job-name=prep_cQTL
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/%x_%j.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=649G

set -euo pipefail
module load BCFtools

[[ $# -ge 1 ]] || { echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch prep_cQTL.sh TISSUE" >&2; exit 1; }
EXPERIMENT="${EXPERIMENT:-target_ALS}"
TISSUE="$1"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')
[[ "$EXPERIMENT" == target_ALS || "$EXPERIMENT" == answer_ALS ]] || { echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2; exit 1; }

BASE="/home/zw529/donglab/data/$EXPERIMENT"
QTL_DIR="$BASE/QTL"
PLINK="$QTL_DIR/plink"
CIRC_MATRIX="$QTL_DIR/circ_matrix.txt"
RAW="$PLINK/joint_all_chrs_matrixEQTL.raw"
BIM="$PLINK/joint_all_chrs_filtered_bed.bim"
PSAM="$PLINK/joint_all_chrs_filtered.psam"
PCA="$PLINK/joint_pca.eigenvec"
MAP_FILE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/ncbi_to_ucsc.txt"
DBSNP_VCF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/GCF_000001405.40.gz"

if [[ "$EXPERIMENT" == target_ALS ]]; then
    OUTDIR="$BASE/$TISSUE_DIR/cQTL"
    METADATA_CSV="$BASE/targetALS_rnaseq_metadata.csv"
    BREAKDOWN_TSV="$QTL_DIR/patient_tissue_breakdown.tsv"
    COV_FILE="$QTL_DIR/covariates.tsv"
    SEX_MISMATCH_FILE="$OUTDIR/potential_sex_mismatches.txt"
else
    [[ "$TISSUE_DIR" == iPSC_Motor_Neuron ]] || { echo "[ERROR] AnswerALS branch supports iPSC_Motor_Neuron only." >&2; exit 1; }
    OUTDIR="$BASE/iPSC_Motor_Neuron/cQTL"
    METADATA_CSV="$BASE/aals_dataportal_datatable.csv"
    BREAKDOWN_TSV=""
    COV_FILE=""
    SEX_MISMATCH_FILE="$BASE/iPSC_Motor_Neuron/eQTL/potential_sex_mismatches.txt"
fi
mkdir -p "$OUTDIR"

for f in "$CIRC_MATRIX" "$RAW" "$BIM" "$PCA" "$MAP_FILE" "$DBSNP_VCF" "$METADATA_CSV"; do
    [[ -s "$f" ]] || { echo "[ERROR] Missing or empty input: $f" >&2; exit 1; }
done
if [[ "$EXPERIMENT" == answer_ALS ]]; then
    [[ -s "$PSAM" ]] || { echo "[ERROR] Missing or empty AnswerALS PSAM: $PSAM" >&2; exit 1; }
fi

echo "============================================"
echo "  $EXPERIMENT cQTL Prep for $TISSUE"
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

echo "Intersecting metadata, genotype, and circRNA matrix..."
python3 << EOF
import gc
import os
import pandas as pd
import numpy as np

experiment, tissue = "$EXPERIMENT", "$TISSUE"

# Keep the original TargetALS import path. For AnswerALS, read the small PSAM roster instead of loading the enormous, very wide PLINK .raw matrix in pandas.
if experiment == "target_ALS":
    raw = pd.read_csv(
        "$RAW",
        sep=r"\\s+",
        usecols=lambda c: c not in ["PAT", "MAT", "SEX", "PHENOTYPE"]
    )
    raw_sample_ids = set(raw["IID"].astype(str))
else:
    psam = pd.read_csv("$PSAM", sep=r"\\s+")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    if "IID" not in psam.columns:
        raise SystemExit("AnswerALS PSAM has no IID column.")
    raw_sample_ids = set(psam["IID"].astype(str))

pca = pd.read_csv("$PCA", sep=r"\\s+").rename(columns={"#IID": "IID"})
headers = pd.read_csv("$CIRC_MATRIX", sep="\\t", nrows=0).columns.tolist()
if "circ_id" not in headers:
    raise SystemExit("Circ matrix must have a circ_id column.")

mismatches = set()
if os.path.isfile("$SEX_MISMATCH_FILE"):
    mismatches = {x.strip() for x in open("$SEX_MISMATCH_FILE") if x.strip()}

coord_to_rsid = {}
for line in open("$TMP_RSID_MAP"):
    f = line.rstrip("\\n").split("\\t")
    if len(f) == 2 and f[1] not in ("", "."):
        coord_to_rsid[f[0].upper()] = f[1]

ucsc_to_ncbi = {}
for line in open("$MAP_FILE"):
    if line.strip():
        ncbi, ucsc = line.strip().split("\\t")
        ucsc_to_ncbi[ucsc.upper()] = ncbi.upper()

if experiment == "target_ALS":
    breakdown = pd.read_csv("$BREAKDOWN_TSV", sep="\\t") if tissue not in ("Cortex", "Spinal_Cord") else None
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
    if tissue in ("Cortex", "Spinal_Cord"):
        mt = meta[meta["tissue"].astype("string").str.contains(pattern, case=False, na=False)].copy()
    else:
        approved = breakdown[breakdown["tissues"].str.contains(tissue, na=False)]["subject_id"].unique()
        mf = meta[meta["externalsubjectid"].isin(approved)].copy()
        mt = mf[mf.apply(lambda r: r.astype(str).str.contains(pattern, case=False).any(), axis=1)].copy()

    def rin(r):
        try:
            return max(
                float(r.iloc[16]) if pd.notnull(r.iloc[16]) else 0,
                float(r.iloc[17]) if pd.notnull(r.iloc[17]) else 0
            )
        except Exception:
            return 0

    mt["RIN_score"] = mt.apply(rin, axis=1)
    mt = mt[(mt["post_mortem_interval_in_hours"] <= 40) & (mt["RIN_score"] >= 3)]
    mu = mt.sort_values("RIN_score", ascending=False).drop_duplicates("externalsubjectid").copy()
    mu["donor_id"] = mu["externalsubjectid"].astype(str)
    mu["matrix_id"] = mu["externalsampleid"].astype(str).str.replace("-", "_", regex=False)
else:
    if tissue != "iPSC_Motor_Neuron":
        raise SystemExit("AnswerALS branch supports iPSC_Motor_Neuron only.")
    meta = pd.read_csv("$METADATA_CSV", low_memory=False, encoding="utf-8-sig")
    id_col = meta.columns[0]
    sex_col = next((c for c in meta.columns if c.strip().lower() == "sex"), None)
    group_col = next((c for c in meta.columns if c.strip().lower() in ("subject_group", "subject group")), None)
    if sex_col is None or group_col is None:
        raise SystemExit("AnswerALS metadata must contain Sex and SUBJECT_GROUP.")
    meta = meta.rename(columns={id_col: "donor_id", sex_col: "sex", group_col: "subject_group"})
    meta["donor_id"] = meta["donor_id"].astype(str).str.strip()
    meta = meta[(meta["donor_id"] != "") & (~meta["donor_id"].isin(mismatches))]
    mu = meta.drop_duplicates("donor_id").copy()
    mu["matrix_id"] = mu["donor_id"]

common = sorted(set(mu["donor_id"].astype(str)) & raw_sample_ids)
final = []
sub_to_matrix = {}
for donor in common:
    prefix = mu.loc[mu["donor_id"].astype(str) == donor, "matrix_id"].iloc[0]
    matches = [c for c in headers if (c == prefix if experiment == "answer_ALS" else c.startswith(prefix))]
    if matches:
        final.append(donor)
        sub_to_matrix[donor] = matches[0]
if not final:
    raise SystemExit("No donors matched metadata, genotype, and circ matrix.")

circ = pd.read_csv(
    "$CIRC_MATRIX",
    sep="\\t",
    usecols=["circ_id"] + [sub_to_matrix[s] for s in final],
    index_col=0
)
circ = circ[[sub_to_matrix[s] for s in final]]
circ.columns = final

def to_rsid(v):
    clean = v.rsplit("_", 1)[0] if "_" in v else v
    p = str(clean).split(":")
    if len(p) < 2:
        return clean
    chrom, pos = p[0], p[1]
    ucsc = chrom if chrom.upper().startswith("CHR") else "chr" + chrom
    ncbi = ucsc_to_ncbi.get(ucsc.upper(), ucsc.upper())
    return coord_to_rsid.get(f"{ncbi}:{pos}".upper(), clean)

snp_out = "$OUTDIR/snp_${TISSUE_DIR}.txt"

if experiment == "target_ALS":
    # Preserve the existing TargetALS SNP-filtering behavior.
    snps = raw[raw["IID"].astype(str).isin(final)].set_index("IID").drop(columns=["FID"]).loc[final].T
    numeric = snps.apply(pd.to_numeric, errors="coerce")
    called = numeric.notna().sum(axis=1)
    af = numeric.sum(axis=1, skipna=True) / (2 * called)
    maf = np.minimum(af, 1 - af)
    snps = snps.loc[(called > 0) & (maf >= 0.05)]

    def groups_pass(row):
        n = pd.to_numeric(row, errors="coerce").dropna().value_counts()
        return len(n) >= 2 and n.min() >= 5

    snps = snps.loc[snps.apply(groups_pass, axis=1)]
    snps.index = [to_rsid(str(x)) for x in snps.index]
    snps = snps.loc[~snps.index.duplicated(keep="first")]
    snps.to_csv(snp_out, sep="\\t", index=True, index_label="snpid", na_rep="NA")
    snps_count = len(snps)

else:
    # AnswerALS-only chunked SNP path:
    # Store hard-call dosages as int8 on disk (-1 means missing), then apply the original MAF and >=5-per-genotype-group rules to bounded variant chunks.
    print("AnswerALS: streaming raw genotypes into a temporary disk-backed matrix...")
    with open("$RAW", "r") as raw_fh:
        raw_header = raw_fh.readline().split()

    if len(raw_header) <= 6 or raw_header[1].lstrip("#") != "IID":
        raise SystemExit("Unexpected PLINK .raw header; expected FID IID PAT MAT SEX PHENOTYPE followed by variants.")

    variant_ids = raw_header[6:]
    n_variants = len(variant_ids)
    n_samples = len(final)
    final_index = {sample: i for i, sample in enumerate(final)}
    temp_matrix = os.path.join("$OUTDIR", f".answerALS_cQTL_genotypes_{os.getpid()}.int8")
    temp_output = snp_out + f".tmp.{os.getpid()}"
    matrix = None
    output_installed = False

    try:
        matrix = np.memmap(temp_matrix, dtype=np.int8, mode="w+", shape=(n_samples, n_variants))
        found = set()

        with open("$RAW", "r") as raw_fh:
            next(raw_fh)  # Skip the wide header.
            for line in raw_fh:
                parts = line.split(None, 6)
                if len(parts) < 7:
                    raise SystemExit("Unexpected short row in PLINK .raw file.")
                sample = parts[1]
                if sample not in final_index:
                    continue

                values = np.fromstring(parts[6].replace("NA", "-1"), sep=" ", dtype=np.float32)
                if values.size != n_variants:
                    raise SystemExit(
                        f"PLINK .raw row for {sample} has {values.size} genotype fields; expected {n_variants}."
                    )

                missing = ~np.isfinite(values) | (values < 0)
                observed = values[~missing]
                if np.any((observed > 2) | (observed != np.floor(observed))):
                    raise SystemExit(
                        f"Non-hard-call dosage found for {sample}; expected AnswerALS non-imputed values 0/1/2."
                    )

                values[missing] = -1
                matrix[final_index[sample], :] = values.astype(np.int8)
                found.add(sample)

        missing_samples = set(final) - found
        if missing_samples:
            raise SystemExit("Samples missing from .raw rows: " + ", ".join(sorted(missing_samples)))

        matrix.flush()
        print(f"AnswerALS: filtering {n_variants:,} variants in 50,000-variant chunks...")
        keep = np.zeros(n_variants, dtype=bool)
        chunk_size = 50000

        for start in range(0, n_variants, chunk_size):
            end = min(start + chunk_size, n_variants)
            block = matrix[:, start:end]

            n0 = np.count_nonzero(block == 0, axis=0)
            n1 = np.count_nonzero(block == 1, axis=0)
            n2 = np.count_nonzero(block == 2, axis=0)
            called = n0 + n1 + n2
            allele_sum = n1 + 2 * n2

            af = np.divide(
                allele_sum,
                2 * called,
                out=np.zeros_like(allele_sum, dtype=np.float64),
                where=called > 0
            )
            maf = np.minimum(af, 1 - af)

            counts = np.vstack((n0, n1, n2))
            categories = np.count_nonzero(counts > 0, axis=0)
            min_group = np.where(counts > 0, counts, np.iinfo(np.int32).max).min(axis=0)

            keep[start:end] = (called > 0) & (maf >= 0.05) & (categories >= 2) & (min_group >= 5)

            del block, n0, n1, n2, called, allele_sum, af, maf, counts, categories, min_group
            if end % 500000 == 0 or end == n_variants:
                print(f"  filtered through {end:,}/{n_variants:,} variants")
            gc.collect()

        candidate_indices = np.flatnonzero(keep)
        print(f"AnswerALS: {len(candidate_indices):,} variants pass MAF and genotype-group filters.")

        seen_ids = set()
        written = 0
        with open(temp_output, "w") as out_fh:
            out_fh.write("snpid\\t" + "\\t".join(final) + "\\n")

            for start in range(0, len(candidate_indices), 10000):
                positions = candidate_indices[start:start + 10000]
                selected_positions = []
                selected_ids = []

                for pos in positions:
                    sid = to_rsid(str(variant_ids[pos]))
                    if sid not in seen_ids:
                        seen_ids.add(sid)
                        selected_positions.append(pos)
                        selected_ids.append(sid)

                if not selected_positions:
                    continue

                values = np.asarray(matrix[:, selected_positions].T, dtype=np.float32)
                values[values < 0] = np.nan
                chunk = pd.DataFrame(values, columns=final)
                chunk.insert(0, "snpid", selected_ids)
                chunk.to_csv(out_fh, sep="\\t", index=False, header=False, na_rep="NA")
                written += len(selected_ids)

                del values, chunk
                gc.collect()

        os.replace(temp_output, snp_out)
        output_installed = True
        snps_count = written
        print(f"AnswerALS: wrote {snps_count:,} filtered SNPs to {snp_out}")

    finally:
        if matrix is not None:
            matrix.flush()
            del matrix
        if os.path.exists(temp_matrix):
            os.remove(temp_matrix)
        if not output_installed and os.path.exists(temp_output):
            os.remove(temp_output)

if experiment == "target_ALS":
    cf = pd.read_csv("$COV_FILE", sep="\\t")
    cov = cf[["externalsampleid"] + [c for c in cf.columns if c not in meta.columns]].drop_duplicates("externalsampleid")
    aligned = mu[mu["donor_id"].isin(final)].copy()
    merged = pd.merge(pd.merge(aligned, cov, on="externalsampleid"), pca, left_on="donor_id", right_on="IID")
    merged = merged.sort_values("RIN_score", ascending=False).drop_duplicates("donor_id")
    merged["is_als"] = merged["subject_group"].apply(lambda x: 1 if "ALS" in str(x) else 0)
    merged["sex_bin"] = merged["sex"].astype(str).str.lower().map({"male": 1, "female": 0})
    covars = merged.set_index("donor_id").reindex(final)[
        ["sex_bin", "age_at_death", "PC1", "PC2", "PC3", "PC4", "PC5", "is_als"]
    ].T
else:
    aligned = mu[mu["donor_id"].isin(final)].copy()
    merged = pd.merge(aligned, pca, left_on="donor_id", right_on="IID", how="inner")
    if len(merged) != len(final):
        raise SystemExit("Some AnswerALS donors lack genotype PCA rows.")
    merged["is_als"] = merged["subject_group"].astype(str).str.contains(
        r"\\bALS\\b", case=False, na=False
    ).astype(int)
    merged["sex_bin"] = merged["sex"].astype(str).str.lower().map({"male": 1, "female": 0})
    covars = merged.set_index("donor_id").reindex(final)[
        ["sex_bin", "PC1", "PC2", "PC3", "PC4", "PC5", "is_als"]
    ].T

circ.to_csv("$OUTDIR/circ_${TISSUE_DIR}.txt", sep="\\t", index=True, index_label="circ_id")
covars.to_csv(
    "$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt",
    sep="\\t", index=True, index_label="id", na_rep="NA"
)
print(f"SUCCESS: Processed {len(final)} donors; {len(circ)} circ features and {snps_count} SNPs.")
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
    key=toupper(ncbi":"pos)
    id=(key in r_map)?r_map[key]:ncbi":"pos
    final_chr=(toupper(ncbi) in inv_map)?inv_map[toupper(ncbi)]:chrom
    if (!seen[id]++) print id,final_chr,pos
}' "$BIM" >> "$OUTDIR/snp_location.txt"

printf "circ_id\tchr\tleft\tright\n" > "$OUTDIR/circ_location.txt"
awk 'BEGIN{OFS="\t"} NR>1 {circ=$1; split(circ,p,":"); split(p[2],xy,"-"); if (!seen[circ]++) print circ,p[1],xy[1],xy[2]}' \
    "$OUTDIR/circ_${TISSUE_DIR}.txt" >> "$OUTDIR/circ_location.txt"

rm -f "$TMP_RSID_MAP"
echo "Location files complete for $TISSUE."
