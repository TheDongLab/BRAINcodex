#!/usr/bin/env bash
#SBATCH --job-name=run_cQTL
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/run_cQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/run_cQTL_%j.err
#SBATCH --time=5-24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

#SBATCH --partition=week

set -euo pipefail
module load R

# ── Arguments ─────────────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch run_cQTL.sh \"Cerebellum\" [interaction] [sub_dir]"
    exit 1
fi

EXPERIMENT="${EXPERIMENT:-target_ALS}"
if [[ "$EXPERIMENT" != "target_ALS" && "$EXPERIMENT" != "answer_ALS" ]]; then
    echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2
    exit 1
fi

TISSUE="$1"
RUN_TYPE="${2:-standard}"   # Default to standard if second argument is empty
SUB_DIR="${3:-}"            # Optional 3rd argument to process subfolder runs natively
if [[ "$RUN_TYPE" != "standard" && "$RUN_TYPE" != "interaction" ]]; then
    echo "[ERROR] RUN_TYPE must be standard or interaction." >&2
    exit 1
fi
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

if [[ "$EXPERIMENT" == "answer_ALS" && "$TISSUE_DIR" != "iPSC_Motor_Neuron" ]]; then
    echo "[ERROR] AnswerALS runner currently supports iPSC_Motor_Neuron only." >&2
    exit 1
fi

echo "========================================================"
if [ -n "$SUB_DIR" ]; then
    echo "  Matrix cQTL run for : $TISSUE ($SUB_DIR) [$RUN_TYPE mode]"
else
    echo "  Matrix cQTL run for tissue : $TISSUE ($RUN_TYPE mode)"
fi
echo "  Experiment: $EXPERIMENT"
echo "  $(date)"
echo "========================================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE="/home/zw529/donglab/data/$EXPERIMENT"
PIPELINE="${QTL_PIPELINE:-/home/zw529/donglab/pipelines/scripts/QTL}"
PLINK="$BASE/QTL/plink"

# Adjust working directory and naming prefixes if a subfolder run is specified
if [ -n "$SUB_DIR" ]; then
    INDIR="$BASE/$TISSUE_DIR/cQTL/$SUB_DIR"
    FILE_PREFIX="${SUB_DIR}"
else
    INDIR="$BASE/$TISSUE_DIR/cQTL"
    FILE_PREFIX="${TISSUE_DIR}"
fi

# Separate output directory based on run type to prevent overwrites
if [ "$RUN_TYPE" == "interaction" ]; then
    OUTDIR="$INDIR/interaction_results"
else
    OUTDIR="$INDIR/results"
fi
mkdir -p "$OUTDIR"

# Match the filenames from your structure dynamically
SNP_FILE="$INDIR/snp_${FILE_PREFIX}.txt"
# circ matrices store fractions (0–1). Models retain fraction-scale effects;
# plotting/reporting scripts convert to percentages exactly once.
CIRC_FILE="$INDIR/circ_${FILE_PREFIX}.txt"
COV_FILE="$INDIR/covariates_${FILE_PREFIX}_encoded.txt"
CIRC_LOC="$INDIR/circ_location.txt"
SNP_LOC="$INDIR/snp_location.txt"
BIM="$PLINK/joint_all_chrs_filtered_bed.bim"

# Output naming
OUTPUT_PREFIX="$OUTDIR/${FILE_PREFIX}_cQTL"
CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
FDR_THRESH=0.05
TOP_N=1000000

# ── Step 0: Pre-flight ────────────────────────────────────────────────
echo "[0] Verifying file existence and alignment..."
for f in "$SNP_FILE" "$CIRC_FILE" "$COV_FILE" "$CIRC_LOC" "$SNP_LOC" "$BIM"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: File missing or empty: $f"
        exit 1
    fi
done

# Alignment Guard
S_N=$(head -n 1 "$SNP_FILE" | awk -F'\t' '{print NF-1}')
C_N_GENES=$(head -n 1 "$CIRC_FILE" | awk -F'\t' '{print NF-1}')
C_N_COVS=$(head -n 1 "$COV_FILE" | awk -F'\t' '{print NF-1}')

echo "Verification: SNPs($S_N), circRNA($C_N_GENES), Covs($C_N_COVS)"

if [[ "$S_N" -ne "$C_N_GENES" ]] || [[ "$S_N" -ne "$C_N_COVS" ]]; then
    echo "FATAL ERROR: Matrices are not aligned! SNP:$S_N, CIRC:$C_N_GENES, COV:$C_N_COVS"
    exit 1
fi

# Exact identity/order check, including duplicated donor IDs.
python3 - "$SNP_FILE" "$CIRC_FILE" "$COV_FILE" <<'PYALIGN'
import sys
headers=[]
for path in sys.argv[1:]:
    with open(path) as f:
        headers.append(f.readline().rstrip("\n").split("\t")[1:])
if not headers[0] or len(set(headers[0])) != len(headers[0]) or not all(h == headers[0] for h in headers):
    raise SystemExit("Donor IDs/order differ between genotype, circRNA, and covariate matrices")
PYALIGN

# ── Step 1: Run Matrix cQTL ───────────────────────────────────────────
echo "[1] Running Matrix cQTL..."
if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_cQTL_LINEAR_CROSS.R" \
        "$SNP_FILE" "$CIRC_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$CIRC_LOC" "$SNP_LOC"
else
    Rscript "$PIPELINE/_cQTL.R" \
        "$SNP_FILE" "$CIRC_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$CIRC_LOC" "$SNP_LOC"
fi

# ── Step 2: Post-processing ──────────────────────────────────────────
echo "[2] Post-processing..."
if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_cQTL_postprocess_LINEAR_CROSS.R" \
        "$CIS_FILE" "$SNP_LOC" "$CIRC_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"
else
    Rscript "$PIPELINE/_cQTL_postprocess.R" \
        "$CIS_FILE" "$SNP_LOC" "$CIRC_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"
fi

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top_for_boxplot.txt"

# [Step 2.5 Omitted: circRNA identifiers do not undergo symbol naming conversion]

# ── Step 3: Manhattan Plot ────────────────────────────────────────────
echo "[3] Generating Manhattan plot..."
if [ "$RUN_TYPE" == "interaction" ]; then
    STD_ANNOTATED="$INDIR/results/${FILE_PREFIX}_cQTL.full_annotated.txt"

    if [ ! -f "$STD_ANNOTATED" ]; then
        echo "FATAL ERROR: Interaction plotting requires baseline standard cQTL results."
        echo "Please run this tissue in standard mode first: sbatch run_cQTL.sh \"$TISSUE\""
        exit 1
    fi

    Rscript "$PIPELINE/_cQTL_manhattan.R" \
        "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH" \
        "$RUN_TYPE" "$STD_ANNOTATED"
else
    Rscript "$PIPELINE/_cQTL_manhattan.R" \
        "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH"
fi

# ── Step 3.5: Regional Locus Zoom ─────────────────────────────────────
if [ "$EXPERIMENT" = "answer_ALS" ]; then
    echo "[3.5] Skipping regional locus zoom: helper still uses TargetALS paths."
elif [ -n "$SUB_DIR" ]; then
    echo "[3.5] Skipping regional locus zoom for stratified subfolder runs."
else
    echo "[3.5] Generating regional locus zoom plots..."
    Rscript "$PIPELINE/_cQTL_regional_zoom.R" "$TISSUE_DIR" "$RUN_TYPE"
fi

# ── Step 4: Boxplots ──────────────────────────────────────────────────
echo "[4] Generating boxplots for all sig. SNPs..."
if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_cQTL_boxplot_LINEAR_CROSS.R" \
        "$TOP_PAIRS" "$SNP_FILE" "$CIRC_FILE" "$COV_FILE" "$SNP_LOC" "$OUTDIR" "${FILE_PREFIX}"
else
    Rscript "$PIPELINE/_cQTL_boxplot.R" \
        "$TOP_PAIRS" "$SNP_FILE" "$CIRC_FILE" "$COV_FILE" "$OUTDIR" "${FILE_PREFIX}"
fi

# ── Step 5: Cleanup Directory Sprawl ──────────────────────────────────
if [ -d "${OUTPUT_PREFIX}" ]; then
    echo "[5] Cleaning up redundant subdirectories..."
    mv "${OUTPUT_PREFIX}"/* "$OUTDIR/" 2>/dev/null || true
    rmdir "${OUTPUT_PREFIX}" 2>/dev/null || true
fi

# ── Final summary ─────────────────────────────────────────────────────
echo "============================================"
echo "  Run complete for : $TISSUE ($RUN_TYPE mode)"
echo "  Experiment       : $EXPERIMENT"
echo "  All outputs saved to: $OUTDIR"
echo "  $(date)"
echo "============================================"
