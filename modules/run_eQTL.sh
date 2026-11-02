#!/bin/bash
#SBATCH --job-name=run_eQTL
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/run_eQTL_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/run_eQTL_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail
module load R

# ── Arguments ─────────────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "ERROR: Missing tissue argument."
    echo "Usage: EXPERIMENT=target_ALS|answer_ALS sbatch run_eQTL.sh \"Cerebellum\" [circ] [standard|interaction] [sub_dir]"
    exit 1
fi

EXPERIMENT="${EXPERIMENT:-target_ALS}"
if [[ "$EXPERIMENT" != "target_ALS" && "$EXPERIMENT" != "answer_ALS" ]]; then
    echo "ERROR: EXPERIMENT must be target_ALS or answer_ALS." >&2
    exit 1
fi

TISSUE="$1"
# circ may appear before or after the existing run-type/subdirectory arguments.
MODE=gene
POSITIONAL=()
for arg in "${@:2}"; do
    if [ "$arg" = circ ]; then MODE=circ; else POSITIONAL+=("$arg"); fi
done
RUN_TYPE="${POSITIONAL[0]:-standard}"
SUB_DIR="${POSITIONAL[1]:-}"
if [ "${#POSITIONAL[@]}" -gt 2 ]; then echo "Too many arguments" >&2; exit 1; fi
if [[ "$RUN_TYPE" != "standard" && "$RUN_TYPE" != "interaction" ]]; then
    echo "ERROR: RUN_TYPE must be standard or interaction." >&2
    exit 1
fi

ANALYSIS=eQTL
if [ "$MODE" = circ ]; then ANALYSIS=circ_eQTL; fi
export EQTL_MODE="$MODE"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

if [[ "$EXPERIMENT" == "answer_ALS" && "$TISSUE_DIR" != "iPSC_Motor_Neuron" ]]; then
    echo "ERROR: AnswerALS runner currently supports iPSC_Motor_Neuron only." >&2
    exit 1
fi

echo "========================================================"
if [ -n "$SUB_DIR" ]; then
    echo "  Matrix eQTL run for : $TISSUE ($SUB_DIR) [$RUN_TYPE mode]"
else
    echo "  Matrix eQTL run for tissue : $TISSUE ($RUN_TYPE mode)"
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
    INDIR="$BASE/$TISSUE_DIR/$ANALYSIS/$SUB_DIR"
    FILE_PREFIX="${SUB_DIR}"
else
    INDIR="$BASE/$TISSUE_DIR/$ANALYSIS"
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
EXPR_FILE="$INDIR/expression_${FILE_PREFIX}.txt"
COV_FILE="$INDIR/covariates_${FILE_PREFIX}_encoded.txt"
GENE_LOC="$INDIR/gene_location.txt"
SNP_LOC="$INDIR/snp_location.txt"
BIM="$PLINK/joint_all_chrs_filtered_bed.bim"

# Output naming (Modified to keep files flat in $OUTDIR)
OUTPUT_PREFIX="$OUTDIR/${FILE_PREFIX}_eQTL"
CIS_FILE="${OUTPUT_PREFIX}.cis.txt"
FDR_THRESH=0.05
TOP_N=1000000

# ── Step 0: Pre-flight ────────────────────────────────────────────────
echo "[0] Verifying file existence and alignment..."
for f in "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$GENE_LOC" "$SNP_LOC" "$BIM"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: File missing or empty: $f"
        exit 1
    fi
done

# Alignment Guard
S_N=$(head -n 1 "$SNP_FILE" | awk -F'\t' '{print NF-1}')
E_N=$(head -n 1 "$EXPR_FILE" | awk -F'\t' '{print NF-1}')
C_N=$(head -n 1 "$COV_FILE" | awk -F'\t' '{print NF-1}')

echo "Verification: SNPs($S_N), Expr($E_N), Covs($C_N)"

if [[ "$S_N" -ne "$E_N" ]] || [[ "$S_N" -ne "$C_N" ]]; then
    echo "FATAL ERROR: Matrices are not aligned! SNP:$S_N, EXPR:$E_N, COV:$C_N"
    exit 1
fi

# Exact identity/order check, including duplicated donor IDs.
python3 - "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" <<'PYALIGN'
import sys
headers=[]
for path in sys.argv[1:]:
    with open(path) as f:
        headers.append(f.readline().rstrip("\n").split("\t")[1:])
if not headers[0] or len(set(headers[0])) != len(headers[0]) or not all(h == headers[0] for h in headers):
    raise SystemExit("Donor IDs/order differ between genotype, phenotype and covariates")
PYALIGN

# ── Step 1: Run Matrix eQTL ───────────────────────────────────────────
echo "[1] Running Matrix eQTL..."
if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_eQTL_LINEAR_CROSS.R" \
        "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$GENE_LOC" "$SNP_LOC"
else
    Rscript "$PIPELINE/_eQTL.R" \
        "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$OUTPUT_PREFIX" "$GENE_LOC" "$SNP_LOC"
fi

# ── Step 2: Post-processing ──────────────────────────────────────────
echo "[2] Post-processing..."
if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_eQTL_postprocess_LINEAR_CROSS.R" \
        "$CIS_FILE" "$SNP_LOC" "$GENE_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"
else
    Rscript "$PIPELINE/_eQTL_postprocess.R" \
        "$CIS_FILE" "$SNP_LOC" "$GENE_LOC" "$OUTPUT_PREFIX" "$FDR_THRESH" "$TOP_N"
fi

ANNOTATED_FILE="${OUTPUT_PREFIX}.full_annotated.txt"
LEAD_FILE="${OUTPUT_PREFIX}.lead_snps.txt"
TOP_PAIRS="${OUTPUT_PREFIX}.top_for_boxplot.txt"

# ── Step 2.5: In-Place Gene Name Conversion (AnnotationHub) ───────────
echo "[2.5] Overwriting Ensembl IDs with common symbols..."
if [ "$MODE" = circ ]; then
    echo "Keeping circRNA IDs and CIRCexplorer2 host annotations."
elif [ -n "$SUB_DIR" ] || [ "$EXPERIMENT" = "answer_ALS" ]; then
    # Use an explicit output directory for stratified and AnswerALS runs.
    TMP_CONVERT="$OUTDIR/tmpscript.R"
    sed '/# TARGET EXECUTION PATHS/,$d' "$PIPELINE/convert_eqtl_names.R" > "$TMP_CONVERT"

    Rscript - "$OUTDIR" "$FILE_PREFIX" "$TMP_CONVERT" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
prefix  <- paste0(args[2], "_eQTL")
source(args[3])

convert_eqtl_genes(file.path(out_dir, paste0(prefix, ".cis.txt")), is_boxplot_file=FALSE)
convert_eqtl_genes(file.path(out_dir, paste0(prefix, ".full_annotated.txt")), is_boxplot_file=FALSE)
convert_eqtl_genes(file.path(out_dir, paste0(prefix, ".FDR0.05.txt")), is_boxplot_file=FALSE)
convert_eqtl_genes(file.path(out_dir, paste0(prefix, ".lead_snps.txt")), is_boxplot_file=FALSE)
convert_eqtl_genes(file.path(out_dir, paste0(prefix, ".top_for_boxplot.txt")), is_boxplot_file=TRUE)
EOF
    rm -f "$TMP_CONVERT"
else
    # Standard TargetALS runs use the existing CLI routine.
    Rscript "$PIPELINE/convert_eqtl_names.R" "$TISSUE_DIR" "$RUN_TYPE"
fi

# ── Step 3: Manhattan Plot ────────────────────────────────────────────
# Guarded: Skip genome-wide Manhattan plots for targeted single-chromosome runs
echo "[3] Generating Manhattan plot..."
if [ "$RUN_TYPE" == "interaction" ]; then
    STD_ANNOTATED="$INDIR/results/${FILE_PREFIX}_eQTL.full_annotated.txt"

    if [ ! -f "$STD_ANNOTATED" ]; then
        echo "FATAL ERROR: Interaction plotting requires baseline standard eQTL results."
        echo "Please run this tissue in standard mode first: sbatch run_eQTL.sh \"$TISSUE\""
        exit 1
    fi

    Rscript "$PIPELINE/_eQTL_manhattan.R" \
        "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH" \
        "$RUN_TYPE" "$STD_ANNOTATED"
else
    Rscript "$PIPELINE/_eQTL_manhattan.R" \
        "$ANNOTATED_FILE" "$LEAD_FILE" "$OUTPUT_PREFIX" "$FDR_THRESH"
fi

# ── Step 3.5: Regional Locus Zoom ─────────────────────────────────────
if [ "$EXPERIMENT" = "answer_ALS" ]; then
    echo "[3.5] Skipping regional locus zoom: helper still uses TargetALS paths."
elif [ -n "$SUB_DIR" ]; then
    echo "[3.5] Skipping standard regional locus zoom for stratified subfolder runs."
else
    echo "[3.5] Generating regional locus zoom plots..."
    Rscript "$PIPELINE/_eQTL_regional_zoom.R" "$TISSUE_DIR" "$RUN_TYPE"
fi

# ── Step 4: Boxplots ──────────────────────────────────────────────────
echo "[4] Generating boxplots for all sig. SNPs..."
if [ "$RUN_TYPE" == "interaction" ]; then
    Rscript "$PIPELINE/_eQTL_boxplot_LINEAR_CROSS.R" \
        "$TOP_PAIRS" "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$SNP_LOC" "$OUTDIR" "${FILE_PREFIX}"
else
    Rscript "$PIPELINE/_eQTL_boxplot.R" \
        "$TOP_PAIRS" "$SNP_FILE" "$EXPR_FILE" "$COV_FILE" "$SNP_LOC" "$OUTDIR" "${FILE_PREFIX}"
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
