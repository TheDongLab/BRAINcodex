#!/bin/bash
#SBATCH --job-name=run_sQTL
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/run_sQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/run_sQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail
module load R

# ── Arguments ─────────────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch run_sQTL.sh \"Cerebellum\" [standard|interaction]"
    exit 1
fi

EXPERIMENT="${EXPERIMENT:-target_ALS}"
if [[ "$EXPERIMENT" != "target_ALS" && "$EXPERIMENT" != "answer_ALS" ]]; then
    echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2
    exit 1
fi

TISSUE="$1"
RUN_TYPE="${2:-standard}"
SUB_DIR="${SUB_DIR:-}"

if [[ "$RUN_TYPE" != "standard" && "$RUN_TYPE" != "interaction" ]]; then
    echo "ERROR: RUN_TYPE must be 'standard' or 'interaction'."
    exit 1
fi

TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

if [[ "$EXPERIMENT" == "answer_ALS" && "$TISSUE_DIR" != "iPSC_Motor_Neuron" ]]; then
    echo "[ERROR] AnswerALS runner currently supports iPSC_Motor_Neuron only." >&2
    exit 1
fi

echo "============================================"
echo "  Matrix sQTL run for tissue : $TISSUE"
echo "  Experiment                : $EXPERIMENT"
echo "  Mode                       : $RUN_TYPE"
echo "  $(date)"
echo "============================================"

# ── Paths ─────────────────────────────────────────────────────────────
BASE="/home/zw529/donglab/data/$EXPERIMENT"
PIPELINE="${QTL_PIPELINE:-/home/zw529/donglab/pipelines/scripts/QTL}"
PLINK="$BASE/QTL/plink"

INDIR="$BASE/$TISSUE_DIR/sQTL"
COV_FILE="$INDIR/covariates_${TISSUE_DIR}_encoded.txt"

if [ "$RUN_TYPE" == "interaction" ]; then
    OUTDIR="$INDIR/interaction_results"
else
    OUTDIR="$INDIR/results"
fi

mkdir -p "$OUTDIR"

# Input matrices
SNP_FILE="$INDIR/snp_${TISSUE_DIR}.txt"
SPLICING_FILE="$INDIR/splicing_${TISSUE_DIR}.txt"
SPLICING_LOC="$INDIR/splicing_location.txt"
SNP_LOC="$INDIR/snp_location.txt"
BIM="$PLINK/joint_all_chrs_filtered_bed.bim"

# Output files
OUTPUT_PREFIX="$OUTDIR/${TISSUE_DIR}_sQTL"
CIS_FILE="${OUTPUT_PREFIX}.cis.txt"

FDR_THRESH=0.05
TOP_N=1000000

# ── Step 0: Pre-flight & Location Harmonization ───────────────────────
echo "[0] Verifying file existence and alignment..."

for f in "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$SPLICING_LOC" "$SNP_LOC" "$BIM"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: File missing or empty: $f"
        exit 1
    fi
done

SPLICING_LOC_BAK="${SPLICING_LOC}.bak"
if [[ ! -f "$SPLICING_LOC_BAK" ]]; then
    cp "$SPLICING_LOC" "$SPLICING_LOC_BAK"
fi

if [[ "$EXPERIMENT" == "answer_ALS" ]]; then
    echo "[0.1] Matching splicing chromosome labels to SNP location labels..."

    SNP_CHR_STYLE=$(awk -F'\t' 'NR>1 && $2!="" {if ($2 ~ /^chr/) print "chr"; else print "plain"; exit}' "$SNP_LOC")
    if [[ -z "$SNP_CHR_STYLE" ]]; then
        echo "[ERROR] Could not determine chromosome naming from $SNP_LOC" >&2
        exit 1
    fi

    awk -F'\t' -v OFS='\t' -v style="$SNP_CHR_STYLE" '
        NR==1 { print; next }
        {
            sub(/^chr/, "", $2)
            if (style == "chr") $2 = "chr" $2
            print
        }
    ' "$SPLICING_LOC_BAK" > "${SPLICING_LOC}.tmp"
    mv "${SPLICING_LOC}.tmp" "$SPLICING_LOC"

    SHARED_CHR=$(comm -12 \
        <(awk -F'\t' 'NR>1 {print $2}' "$SNP_LOC" | sort -u) \
        <(awk -F'\t' 'NR>1 {print $2}' "$SPLICING_LOC" | sort -u) | wc -l)

    echo "      -> Shared chromosome labels: $SHARED_CHR"
    if [[ "$SHARED_CHR" -eq 0 ]]; then
        echo "[ERROR] SNP and splicing location files have no matching chromosome labels." >&2
        exit 1
    fi
else
    # Preserve the existing TargetALS location-label cleanup.
    echo "[0.1] Standardizing chromosome labels in splicing location file..."
    sed -E 's/\t(chr[0-9XY]+):([+-])\t/\t\1\t/g' \
        "$SPLICING_LOC_BAK" > "$SPLICING_LOC"
    echo "      -> Chromosome strings cleaned in $SPLICING_LOC"
fi

# Alignment guard
S_N=$(head -n 1 "$SNP_FILE" | awk -F'\t' '{print NF-1}')
SPL_N=$(head -n 1 "$SPLICING_FILE" | awk -F'\t' '{print NF-1}')
C_N=$(head -n 1 "$COV_FILE" | awk -F'\t' '{print NF-1}')

echo "Verification: SNPs($S_N), Splicing($SPL_N), Covs($C_N)"

if [[ "$S_N" -ne "$SPL_N" ]] || [[ "$S_N" -ne "$C_N" ]]; then
    echo "FATAL ERROR: Matrices are not aligned! SNP:$S_N, SPLICING:$SPL_N, COV:$C_N"
    exit 1
fi

python3 - "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" <<'PYALIGN'
import sys

headers = []
for path in sys.argv[1:]:
    with open(path) as f:
        headers.append(f.readline().rstrip("\n").split("\t")[1:])

if (
    not headers[0]
    or len(set(headers[0])) != len(headers[0])
    or not all(header == headers[0] for header in headers)
):
    raise SystemExit(
        "Donor IDs/order differ between genotype, splicing, and covariate matrices"
    )
PYALIGN

# ── Step 1: Run Matrix sQTL ───────────────────────────────────────────
echo "[1] Running Matrix sQTL..."

if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_sQTL_LINEAR_CROSS.R" \
        "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$SPLICING_LOC" "$SNP_LOC"
else
    Rscript "$PIPELINE/_sQTL.R" \
        "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$SPLICING_LOC" "$SNP_LOC"
fi

# ── Step 2: Post-processing ───────────────────────────────────────────
echo "[2] Post-processing..."

if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_sQTL_postprocess_LINEAR_CROSS.R" \
        "$CIS_FILE" "$SNP_LOC" "$SPLICING_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"
else
    Rscript "$PIPELINE/_sQTL_postprocess.R" \
        "$CIS_FILE" "$SNP_LOC" "$SPLICING_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"
fi

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top_for_boxplot.txt"

# ── Step 3: Manhattan Plot ────────────────────────────────────────────
echo "[3] Generating Manhattan plot..."

if [ "$RUN_TYPE" == "interaction" ]; then
    STD_ANNOTATED="$INDIR/results/${TISSUE_DIR}_sQTL.full_annotated.txt"

    if [ ! -f "$STD_ANNOTATED" ]; then
        echo "FATAL ERROR: Interaction plotting requires baseline standard sQTL results."
        echo "Run standard mode first:"
        echo "sbatch $0 \"$TISSUE\" standard"
        exit 1
    fi

    Rscript "$PIPELINE/_sQTL_manhattan.R" \
        "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH" "$RUN_TYPE" "$STD_ANNOTATED"
else
    Rscript "$PIPELINE/_sQTL_manhattan.R" \
        "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH"
fi

# ── Step 3.5: Regional Locus Zoom ─────────────────────────────────────
if [ "$EXPERIMENT" = "answer_ALS" ]; then
    echo "[3.5] Generating AnswerALS regional locus zoom plots..."
    Rscript "$PIPELINE/_sQTL_regional_zoom.R" "$TISSUE_DIR" "$RUN_TYPE"
elif [ -n "$SUB_DIR" ]; then
    echo "[3.5] Skipping standard regional locus zoom for stratified subfolder runs."
else
    echo "[3.5] Generating regional locus zoom plots..."
    Rscript "$PIPELINE/_sQTL_regional_zoom.R" "$TISSUE_DIR" "$RUN_TYPE"
fi

# ── Step 4: Boxplots ──────────────────────────────────────────────────
echo "[4] Generating boxplots for all significant SNPs..."

if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_sQTL_boxplot_LINEAR_CROSS.R" \
        "$TOP_PAIRS" "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$SNP_LOC" "$OUTDIR" "$TISSUE_DIR"
else
    Rscript "$PIPELINE/_sQTL_boxplot.R" \
        "$TOP_PAIRS" "$SNP_FILE" "$SPLICING_FILE" "$COV_FILE" "$SNP_LOC" "$OUTDIR" "$TISSUE_DIR"
fi

# ── Step 5: Cleanup Directory Sprawl ──────────────────────────────────
if [ -d "${OUTPUT_PREFIX}" ]; then
    echo "[5] Cleaning up redundant subdirectories..."
    mv "${OUTPUT_PREFIX}"/* "$OUTDIR/" 2>/dev/null || true
    rmdir "${OUTPUT_PREFIX}" 2>/dev/null || true
fi

# ── Final Summary ─────────────────────────────────────────────────────
echo "============================================"
echo "  Run complete for : $TISSUE (Mode: $RUN_TYPE)"
echo "  Experiment       : $EXPERIMENT"
echo "  All outputs saved to: $OUTDIR"
echo "  $(date)"
echo "============================================"
