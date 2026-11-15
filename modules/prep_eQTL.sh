#!/bin/bash
#SBATCH --job-name=prep_eQTL
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=649G

set -euo pipefail

module purge
unset PYTHONPATH
module load Python/3.12.3-GCCcore-13.3.0
module load Python-bundle-PyPI/2024.06-GCCcore-13.3.0
module load BCFtools

if [ "$#" -lt 1 ]; then
    echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch prep_eQTL.sh TISSUE [circ]" >&2
    exit 1
fi

EXPERIMENT="${EXPERIMENT:-target_ALS}"
TISSUE="$1"
MODE="${2:-gene}"

case "$EXPERIMENT" in
    target_ALS|answer_ALS) ;;
    *) echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2; exit 1 ;;
esac

case "$MODE" in
    gene|circ) ;;
    *) echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch prep_eQTL.sh TISSUE [circ]" >&2; exit 1 ;;
esac

ANALYSIS=eQTL
if [ "$MODE" = circ ]; then ANALYSIS=circ_eQTL; fi

SCRIPT_DIR="${QTL_PIPELINE:-/home/zw529/donglab/pipelines/scripts/QTL}"
if [ "$MODE" = circ ]; then
    module load Python/3.12.3-GCCcore-13.3.0 Python-bundle-PyPI/2024.06-GCCcore-13.3.0
fi

TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')
BASE="/home/zw529/donglab/data/$EXPERIMENT"

if [ "$EXPERIMENT" = "target_ALS" ]; then
    OUTDIR="$BASE/$TISSUE_DIR/$ANALYSIS"
    QTL_DIR="$BASE/QTL"
    METADATA_CSV="$BASE/targetALS_rnaseq_metadata.csv"
    BREAKDOWN_TSV="$QTL_DIR/patient_tissue_breakdown.tsv"
    COV_FILE="$QTL_DIR/covariates.tsv"
else
    if [ "$TISSUE_DIR" != "iPSC_Motor_Neuron" ]; then
        echo "[ERROR] This AnswerALS branch currently supports iPSC_Motor_Neuron only." >&2
        exit 1
    fi
    OUTDIR="$BASE/iPSC_Motor_Neuron/$ANALYSIS"
    QTL_DIR="$BASE/QTL"
    METADATA_CSV="$BASE/aals_dataportal_datatable.csv"
    BREAKDOWN_TSV=""
    COV_FILE=""
fi

mkdir -p "$OUTDIR"

### Paths ###
PLINK="$QTL_DIR/plink"
EXPR="$QTL_DIR/expression_matrix.txt"
if [ "$MODE" = circ ]; then EXPR="$QTL_DIR/circ_rpm_matrix.txt"; fi
SEX_MISMATCH_FILE="$OUTDIR/potential_sex_mismatches.txt"

# --- GENOTYPE INPUTS ---
RAW="$PLINK/joint_all_chrs_matrixEQTL.raw"
BIM="$PLINK/joint_all_chrs_filtered_bed.bim"
PSAM="$PLINK/joint_all_chrs_filtered.psam"
PCA="$PLINK/joint_pca.eigenvec"

# --- REFERENCE FILES ---
GENCODE_DIR="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode"
GTF_BED6="$GENCODE_DIR/gencode.v49.annotation.gene.bed6"
MAP_FILE="$GENCODE_DIR/ncbi_to_ucsc.txt"
DBSNP_VCF="$GENCODE_DIR/GCF_000001405.40.gz"

for required in "$METADATA_CSV" "$EXPR" "$RAW" "$BIM" "$PCA" "$GTF_BED6" "$MAP_FILE" "$DBSNP_VCF"; do
    [[ -s "$required" ]] || {
        echo "[ERROR] Required input is missing or empty: $required" >&2
        exit 1
    }
done
if [ "$EXPERIMENT" = "answer_ALS" ]; then
    [[ -s "$PSAM" ]] || { echo "[ERROR] Missing or empty AnswerALS PSAM: $PSAM" >&2; exit 1; }
fi

echo "============================================"
echo "  $EXPERIMENT $ANALYSIS Prep for $TISSUE"
echo "============================================"

# ─────────────────────────────────────────────────────────────────
# STEP A: HIGH-PERFORMANCE RSID EXTRACTION VIA BCFTOOLS -T
# ─────────────────────────────────────────────────────────────────
echo "Generating target position filter list..."
TMP_TARGETS="$OUTDIR/tmp_targets_ncbi.txt"

awk -v map_file="$MAP_FILE" '
BEGIN {
    OFS="\t"
    while ((getline < map_file) > 0) { map[toupper($2)] = toupper($1) }
    close(map_file)
}
{
    ucsc = ($1 ~ /^[Cc][Hh][Rr]/) ? toupper($1) : "CHR"toupper($1);
    ncbi = (ucsc in map) ? map[ucsc] : ucsc;
    print ncbi, $4, $4
}' "$BIM" > "$TMP_TARGETS"

echo "Streaming dbSNP variant IDs matching your dataset..."
TMP_RSID_MAP="$OUTDIR/tmp_coord_to_rsid.txt"
bcftools query -T "$TMP_TARGETS" -f '%CHROM:%POS\t%ID\n' "$DBSNP_VCF" > "$TMP_RSID_MAP"
rm -f "$TMP_TARGETS"

# ─────────────────────────────────────────────────────────────────
# STEP B: METADATA INTERSECTION & MATRIX GENERATION
# ─────────────────────────────────────────────────────────────────
python3 << EOF
import gc
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

circ_mode = "$MODE" == "circ"
experiment = "$EXPERIMENT"
tissue = "$TISSUE"
sys.path.insert(0, "$SCRIPT_DIR")

# Preserve the existing TargetALS raw-matrix import. AnswerALS reads the small
# PSAM roster and streams genotype data only when writing the SNP matrix.
if experiment == "target_ALS":
    raw_df = pd.read_csv(
        "$RAW",
        sep=r'\\s+',
        usecols=lambda x: x not in ['PAT', 'MAT', 'SEX', 'PHENOTYPE']
    )
    raw_sample_ids = set(raw_df['IID'].astype(str))
else:
    psam = pd.read_csv("$PSAM", sep=r'\\s+')
    psam.columns = [c.lstrip("#") for c in psam.columns]
    if "IID" not in psam.columns:
        raise SystemExit("AnswerALS PSAM has no IID column.")
    raw_sample_ids = set(psam["IID"].astype(str))

pca = pd.read_csv("$PCA", sep=r'\\s+').rename(columns={'#IID': 'IID'})
expr_headers = pd.read_csv("$EXPR", sep='\\t', nrows=0).columns.tolist()

ucsc_to_ncbi = {}
with open("$MAP_FILE", 'r') as f:
    for line in f:
        if line.strip():
            ncbi, ucsc = line.strip().split('\\t')
            ucsc_to_ncbi[ucsc.upper()] = ncbi.upper()

coord_to_rsid = {}
with open("$TMP_RSID_MAP", 'r') as f:
    for line in f:
        if line.strip():
            coord, rsid = line.strip().split('\\t')
            if rsid and rsid != '.':
                coord_to_rsid[coord.upper().strip()] = rsid.strip()

try:
    with open("$SEX_MISMATCH_FILE") as f:
        sex_mismatch_samples = {line.strip() for line in f if line.strip()}
except FileNotFoundError:
    sex_mismatch_samples = set()

if circ_mode:
    mismatch_path = Path("$BASE/$TISSUE_DIR/eQTL/potential_sex_mismatches.txt")
    if mismatch_path.exists():
        sex_mismatch_samples.update(
            x.strip() for x in mismatch_path.read_text().splitlines() if x.strip()
        )

if experiment == "target_ALS":
    # Keep existing TargetALS metadata and tissue-selection behavior.
    breakdown = None
    if not circ_mode and tissue not in ("Cortex", "Spinal_Cord"):
        breakdown = pd.read_csv("$BREAKDOWN_TSV", sep='\\t')

    meta = pd.read_csv("$METADATA_CSV", low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower()
    meta = meta[~meta['externalsampleid'].isin(sex_mismatch_samples)].copy()

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
    if circ_mode:
        patterns["Lumbar_Spinal_Cord"] += "|Lumbosacaral"
    pattern = patterns.get(tissue, tissue)

    if circ_mode or tissue in ("Cortex", "Spinal_Cord"):
        meta_tissue = meta[meta["tissue"].astype("string").str.contains(pattern, case=False, na=False)].copy()
    else:
        approved_subjects = breakdown[
            breakdown['tissues'].str.contains(tissue, na=False)
        ]['subject_id'].unique()
        meta_filt = meta[meta['externalsubjectid'].isin(approved_subjects)].copy()
        mask = meta_filt.apply(
            lambda row: row.astype(str).str.contains(pattern, case=False).any(), axis=1
        )
        meta_tissue = meta_filt[mask].copy()

    def get_rin(row):
        try:
            val1 = float(row.iloc[16]) if pd.notnull(row.iloc[16]) else 0
            val2 = float(row.iloc[17]) if pd.notnull(row.iloc[17]) else 0
            return max(val1, val2)
        except Exception:
            return 0

    meta_tissue['RIN_score'] = meta_tissue.apply(get_rin, axis=1)
    meta_tissue = meta_tissue[
        (meta_tissue['post_mortem_interval_in_hours'] <= 40) &
        (meta_tissue['RIN_score'] >= 3)
    ]
    meta_unique = meta_tissue.sort_values(
        'RIN_score', ascending=False
    ).drop_duplicates(subset='externalsubjectid')
    subject_col = 'externalsubjectid'
    sample_col = 'externalsampleid'

else:
    if tissue != "iPSC_Motor_Neuron":
        raise SystemExit("AnswerALS branch supports iPSC_Motor_Neuron only.")
    meta = pd.read_csv("$METADATA_CSV", low_memory=False, encoding='utf-8-sig')
    id_col = meta.columns[0]
    sex_col = next((c for c in meta.columns if c.strip().lower() == 'sex'), None)
    group_col = next(
        (c for c in meta.columns if c.strip().lower() in ('subject_group', 'subject group')),
        None
    )
    if sex_col is None:
        raise SystemExit("AnswerALS metadata has no Sex column.")
    if group_col is None:
        raise SystemExit("AnswerALS metadata has no SUBJECT_GROUP column.")

    meta = meta.rename(columns={
        id_col: 'externalsubjectid',
        sex_col: 'sex',
        group_col: 'subject_group'
    })
    meta['externalsubjectid'] = meta['externalsubjectid'].astype(str).str.strip()
    meta['externalsampleid'] = meta['externalsubjectid']
    meta = meta[meta['externalsubjectid'].ne('')].copy()
    meta = meta[~meta['externalsubjectid'].isin(sex_mismatch_samples)].copy()
    meta_unique = meta.drop_duplicates(subset='externalsubjectid').copy()
    meta_unique['RIN_score'] = 0
    subject_col = 'externalsubjectid'
    sample_col = 'externalsampleid'

# Intersect metadata, genotype samples, and expression columns.
if experiment == "target_ALS":
    meta_unique['hra_clean'] = meta_unique['externalsampleid'].astype(str).str.replace('-', '_', regex=False)
else:
    meta_unique['hra_clean'] = meta_unique['externalsubjectid'].astype(str)

common_subjects = sorted(
    set(meta_unique[subject_col].astype(str)) & raw_sample_ids
)

sub_to_hra = {}
final_aligned_subjects = []
for subject in common_subjects:
    row = meta_unique.loc[meta_unique[subject_col].astype(str) == subject].iloc[0]
    hra = row['hra_clean']
    if circ_mode or experiment == "answer_ALS":
        matches = [c for c in expr_headers if c == hra]
    else:
        matches = [c for c in expr_headers if c.startswith(hra)]
    if matches:
        final_aligned_subjects.append(subject)
        sub_to_hra[subject] = matches[0]

if not final_aligned_subjects:
    raise SystemExit("No shared donors remain across metadata, genotype, and expression inputs.")

expr_id = 'circ_id' if circ_mode else 'gene_id'
expr = pd.read_csv(
    "$EXPR", sep='\\t',
    usecols=[expr_id] + [sub_to_hra[s] for s in final_aligned_subjects],
    index_col=0
)
expr_final = expr[[sub_to_hra[s] for s in final_aligned_subjects]]
expr_final.columns = final_aligned_subjects

# Keep the original TargetALS genotype DataFrame path. AnswerALS avoids
# materializing the full variants-by-donors DataFrame.
if experiment == "target_ALS":
    snp_matrix = raw_df[raw_df['IID'].astype(str).isin(final_aligned_subjects)]
    snp_matrix = snp_matrix.set_index('IID').drop(columns=['FID']).loc[final_aligned_subjects]
    snp_final = snp_matrix.T
else:
    snp_final = None

def convert_to_rsid(full_id):
    clean_id = full_id.rsplit('_', 1)[0] if '_' in full_id else full_id
    parts = str(clean_id).split(':')
    if len(parts) >= 2:
        chrom, pos = parts[0], parts[1]
        ucsc_chrom = chrom if chrom.upper().startswith('CHR') else f"chr{chrom}"
        ncbi_chrom = ucsc_to_ncbi.get(ucsc_chrom.upper(), ucsc_chrom.upper())
        coord_key = f"{ncbi_chrom}:{pos}".upper().strip()
        return coord_to_rsid.get(coord_key, clean_id)
    return clean_id

# Preserve the existing TargetALS gene-mode SNP filtering.
if experiment == "target_ALS" and not circ_mode:
    min_maf = 0.05
    snp_numeric = snp_final.apply(pd.to_numeric, errors='coerce')
    n_called = snp_numeric.notna().sum(axis=1)
    allele_freq = snp_numeric.sum(axis=1, skipna=True) / (2 * n_called)
    maf = np.minimum(allele_freq, 1 - allele_freq)
    n_before_maf = len(snp_final)
    snp_final = snp_final.loc[(n_called > 0) & (maf >= min_maf)]
    print(f"DEBUG: Tissue-specific MAF >= {min_maf}: {n_before_maf} -> {len(snp_final)} SNPs.")

    min_genotype_n = 5
    def genotype_groups_pass(row):
        counts = pd.to_numeric(row, errors='coerce').dropna().value_counts()
        return len(counts) >= 2 and counts.min() >= min_genotype_n

    n_before_genotype = len(snp_final)
    snp_final = snp_final.loc[snp_final.apply(genotype_groups_pass, axis=1)]
    print(
        f"DEBUG: Minimum genotype-group N >= {min_genotype_n}: "
        f"{n_before_genotype} -> {len(snp_final)} SNPs."
    )

if experiment == "target_ALS":
    snp_final.index = [convert_to_rsid(str(idx)) for idx in snp_final.index]
    snp_final = snp_final[~snp_final.index.duplicated(keep='first')]

# Covariates.
if experiment == "target_ALS":
    cov_full = pd.read_csv("$COV_FILE", sep='\\t')
    cov_cols_to_keep = ['externalsampleid'] + [
        c for c in cov_full.columns if c not in meta.columns
    ]
    cov = cov_full[cov_cols_to_keep].drop_duplicates(subset='externalsampleid')

    meta_aligned = meta_unique[
        meta_unique['externalsubjectid'].isin(final_aligned_subjects)
    ].copy()
    cov_merged = pd.merge(meta_aligned, cov, on='externalsampleid', how='inner')
    cov_merged = pd.merge(cov_merged, pca, left_on='externalsubjectid', right_on='IID', how='inner')
    cov_merged = cov_merged.sort_values('RIN_score', ascending=False).drop_duplicates(
        subset='externalsubjectid'
    )
    cov_merged['is_als'] = cov_merged['subject_group'].apply(
        lambda x: 1 if 'ALS' in str(x) else 0
    )
    cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map(
        {'male': 1, 'female': 0}
    )
    cov_final = cov_merged.set_index('externalsubjectid').reindex(final_aligned_subjects)
    cov_final = cov_final[
        ['sex_bin', 'age_at_death', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'is_als']
    ].T
else:
    meta_aligned = meta_unique[
        meta_unique['externalsubjectid'].isin(final_aligned_subjects)
    ].copy()
    cov_merged = pd.merge(
        meta_aligned, pca, left_on='externalsubjectid', right_on='IID', how='inner'
    )
    cov_merged['is_als'] = cov_merged['subject_group'].astype(str).str.contains(
        r'\\bALS\\b', case=False, na=False
    ).astype(int)
    cov_merged['sex_bin'] = cov_merged['sex'].astype(str).str.lower().map(
        {'male': 1, 'female': 0}
    )
    cov_final = cov_merged.set_index('externalsubjectid').reindex(final_aligned_subjects)
    cov_final = cov_final[
        ['sex_bin', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'is_als']
    ].T

# Circ mode: TargetALS continues to use the shared helper. For AnswerALS,
# reproduce its circ/sample behavior without constructing a full SNP DataFrame.
if circ_mode and experiment == "target_ALS":
    from circ_eqtl import prepare_circ
    expr_final, snp_final, cov_final = prepare_circ(
        expr_final, snp_final, cov_final, "$OUTDIR", "$TISSUE_DIR",
        "$QTL_DIR/circ_rpm_sample_metadata.tsv", sub_to_hra
    )
    num_final = expr_final.shape[1]

elif circ_mode and experiment == "answer_ALS":
    if expr_final.index.has_duplicates or expr_final.columns.has_duplicates:
        raise ValueError("Duplicated circRNA or donor IDs")

    cov_numeric = cov_final.apply(pd.to_numeric, errors="coerce")
    donors = [
        d for d in expr_final.columns
        if d in cov_numeric.columns and np.isfinite(cov_numeric[d]).all()
    ]
    if len(donors) < 3:
        raise ValueError("Fewer than three donors with complete covariates")

    rpm = expr_final.loc[:, donors]
    cov_final = cov_numeric.loc[:, donors]

    values = rpm.to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("RPM values must be finite and nonnegative")

    variable = rpm.max(axis=1) > rpm.min(axis=1)
    print(
        f"Circ mode: {len(donors)} complete donors; "
        f"{int(variable.sum())}/{len(variable)} variable circRNAs; "
        "no read/prevalence filter"
    )
    rpm = rpm.loc[variable].copy()
    if rpm.empty:
        raise ValueError("No variable circRNAs in the final donor cohort")

    expr_final = pd.DataFrame(
        norm.ppf((rpm.rank(axis=1, method="average") - 0.5) / len(donors)),
        index=rpm.index,
        columns=donors
    )

    coords = rpm.index.to_series().str.extract(
        r'^(?P<chr>[^:]+):(?P<start>\\d+)-(?P<end>\\d+):(?P<strand>[+-])$'
    )
    if coords.isna().any().any():
        raise ValueError("Invalid circRNA coordinate ID")
    left = coords.start.astype("int64") + 1
    right = coords.end.astype("int64")
    if (left > right).any():
        raise ValueError("Invalid circRNA interval")

    loc = pd.DataFrame({
        "geneid": rpm.index,
        "chr": coords.chr.to_numpy(),
        "left": left.to_numpy(),
        "right": right.to_numpy()
    })
    loc.to_csv(Path("$OUTDIR") / "gene_location.txt", sep="\\t", index=False)

    provenance_file = Path("$QTL_DIR/circ_rpm_sample_metadata.tsv")
    provenance = pd.read_csv(provenance_file, sep="\\t", index_col=0)
    genes, transcripts = {}, {}
    selected = set(rpm.index)

    for donor in donors:
        sample = sub_to_hra[donor]
        path = provenance.loc[sample, "circ_file"]
        tab = pd.read_csv(
            path,
            sep="\\t",
            usecols=["chrom", "start", "end", "strand", "geneName", "isoformName"]
        )
        ids = (
            tab.chrom.astype(str) + ":" +
            tab.start.astype(str) + "-" +
            tab.end.astype(str) + ":" +
            tab.strand.astype(str)
        )
        for cid, gene, transcript in zip(ids, tab.geneName, tab.isoformName):
            if cid not in selected:
                continue
            for value, target in [(gene, genes), (transcript, transcripts)]:
                if pd.notna(value) and str(value) not in ("", ".", "NA"):
                    target.setdefault(cid, set()).add(str(value))

    pd.DataFrame({
        "geneid": rpm.index,
        "host_gene": [";".join(sorted(genes.get(c, []))) for c in rpm.index],
        "host_transcript": [";".join(sorted(transcripts.get(c, []))) for c in rpm.index]
    }).to_csv(Path("$OUTDIR") / "circ_annotation.tsv", sep="\\t", index=False)

    rpm.to_csv(
        Path("$OUTDIR") / f"circ_rpm_$TISSUE_DIR.txt",
        sep="\\t",
        index_label="circ_id"
    )
    pd.DataFrame({
        "donor_id": donors,
        "sample_id": [sub_to_hra[d] for d in donors]
    }).to_csv(Path("$OUTDIR") / "selected_samples.tsv", sep="\\t", index=False)

    final_aligned_subjects = donors
    num_final = len(donors)

else:
    num_final = len(final_aligned_subjects)

# AnswerALS SNPs are processed in chunks after circ mode has finalized its donor set.
if experiment == "answer_ALS":
    print("AnswerALS: streaming raw genotypes into a temporary disk-backed matrix...")
    with open("$RAW", "r") as raw_fh:
        raw_header = raw_fh.readline().split()
    if len(raw_header) <= 6 or raw_header[1].lstrip("#") != "IID":
        raise SystemExit("Unexpected PLINK .raw header; expected FID IID PAT MAT SEX PHENOTYPE followed by variants.")

    variant_ids = raw_header[6:]
    n_variants = len(variant_ids)
    n_samples = len(final_aligned_subjects)
    final_index = {sample: i for i, sample in enumerate(final_aligned_subjects)}
    temp_matrix = os.path.join("$OUTDIR", f".answerALS_eQTL_genotypes_{os.getpid()}.int8")
    snp_out = "$OUTDIR/snp_${TISSUE_DIR}.txt"
    temp_output = snp_out + f".tmp.{os.getpid()}"
    matrix = None
    output_installed = False

    try:
        matrix = np.memmap(temp_matrix, dtype=np.int8, mode="w+", shape=(n_samples, n_variants))
        found = set()
        with open("$RAW", "r") as raw_fh:
            next(raw_fh)
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

        missing_samples = set(final_aligned_subjects) - found
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
            out_fh.write("snpid\\t" + "\\t".join(final_aligned_subjects) + "\\n")
            for start in range(0, len(candidate_indices), 10000):
                positions = candidate_indices[start:start + 10000]
                selected_positions, selected_ids = [], []
                for pos in positions:
                    full_id = str(variant_ids[pos])
                    clean_id = full_id.rsplit("_", 1)[0] if "_" in full_id else full_id
                    parts = clean_id.split(":")
                    if len(parts) >= 2:
                        chrom, pos_bp = parts[0], parts[1]
                        ucsc_chrom = chrom if chrom.upper().startswith("CHR") else f"chr{chrom}"
                        ncbi_chrom = ucsc_to_ncbi.get(ucsc_chrom.upper(), ucsc_chrom.upper())
                        coord_key = f"{ncbi_chrom}:{pos_bp}".upper().strip()
                        sid = coord_to_rsid.get(coord_key, clean_id)
                    else:
                        sid = clean_id
                    if sid not in seen_ids:
                        seen_ids.add(sid)
                        selected_positions.append(pos)
                        selected_ids.append(sid)

                if not selected_positions:
                    continue

                values = np.asarray(matrix[:, selected_positions].T, dtype=np.float32)
                values[values < 0] = np.nan
                chunk = pd.DataFrame(values, columns=final_aligned_subjects)
                chunk.insert(0, "snpid", selected_ids)
                chunk.to_csv(out_fh, sep="\\t", index=False, header=False, na_rep="NA")
                written += len(selected_ids)
                del values, chunk
                gc.collect()

        os.replace(temp_output, snp_out)
        output_installed = True
        snp_count = written
        print(f"AnswerALS: wrote {snp_count:,} filtered SNPs to {snp_out}")

    finally:
        if matrix is not None:
            matrix.flush()
            del matrix
        if os.path.exists(temp_matrix):
            os.remove(temp_matrix)
        if not output_installed and os.path.exists(temp_output):
            os.remove(temp_output)

    # Preserve the expression/covariate filenames and circ/gene labels used downstream.
    expr_final.to_csv(
        "$OUTDIR/expression_${TISSUE_DIR}.txt",
        sep="\\t", index=True, index_label="geneid"
    )
    cov_final.to_csv(
        "$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt",
        sep="\\t", index=True, index_label="id", quoting=0, na_rep="NA"
    )
else:
    # Original TargetALS output path.
    snp_final.to_csv(
        "$OUTDIR/snp_${TISSUE_DIR}.txt",
        sep="\\t", index=True, index_label="snpid", na_rep="NA"
    )
    expr_final.to_csv(
        "$OUTDIR/expression_${TISSUE_DIR}.txt",
        sep="\\t", index=True, index_label="geneid"
    )
    cov_final.to_csv(
        "$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt",
        sep="\\t", index=True, index_label="id", quoting=0, na_rep="NA"
    )
    snp_count = len(snp_final)

print(f"SUCCESS: Processed {num_final} unique samples for {tissue}; {snp_count} SNPs.")
EOF

# ─────────────────────────────────────────────────────────────────
# STEP C: GENERATE ALIGNED LOCATION FILES
# ─────────────────────────────────────────────────────────────────
echo "Generating order-preserved location files for $TISSUE..."

printf "snpid\tchr\tpos\n" > "$OUTDIR/snp_location.txt"
awk -v ncbi_map="$MAP_FILE" -v rsid_map="$TMP_RSID_MAP" '
BEGIN {
    OFS="\t"
    while ((getline < ncbi_map) > 0) { n_map[toupper($2)] = toupper($1); inv_map[toupper($1)] = $2 }
    close(ncbi_map)
    while ((getline < rsid_map) > 0) { r_map[toupper($1)] = $2 }
    close(rsid_map)
}
{
    split($2, parts, ":")
    chrom = parts[1]
    pos = parts[2]
    ucsc = (chrom ~ /^[Cc][Hh][Rr]/) ? toupper(chrom) : "CHR"toupper(chrom)
    ncbi = (ucsc in n_map) ? n_map[ucsc] : ucsc
    coord_key = ncbi":"pos
    final_id = (toupper(coord_key) in r_map) ? r_map[toupper(coord_key)] : coord_key
    final_chr = (toupper(ncbi) in inv_map) ? inv_map[toupper(ncbi)] : chrom
    if (!seen[final_id]++) print final_id, final_chr, pos
}' "$BIM" >> "$OUTDIR/snp_location.txt"

if [ "$MODE" != circ ]; then
    printf "geneid\tchr\tleft\tright\n" > "$OUTDIR/gene_location.txt"
    awk 'BEGIN{OFS="\t"} {gene=$4; sub(/___.*$/, "", gene); if (!seen[gene]++) print gene, $1, $2, $3}' \
        "$GTF_BED6" >> "$OUTDIR/gene_location.txt"
fi

rm -f "$TMP_RSID_MAP"
echo "Location files complete for $TISSUE."
