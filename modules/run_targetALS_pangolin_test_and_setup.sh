#!/bin/bash
#SBATCH --job-name=pangolin_targetALS
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/prep_test_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/prep_test_%j.err

set -euo pipefail

PYTHON="$HOME/donglab/pipelines/modules/miniconda3/bin/python"
PANGOLIN="$HOME/donglab/pipelines/modules/miniconda3/bin/pangolin"
SOURCE="$HOME/donglab/pipelines/modules/Pangolin"
OUT="$HOME/donglab/data/target_ALS/QTL/pangolin"
FASTA="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
DB="$SOURCE/references/gencode.v49.all_transcripts.db"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$OUT"

for file in \
    "$PYTHON" \
    "$PANGOLIN" \
    "$SOURCE/scripts/prepare_targetals_variants.py" \
    "$FASTA" \
    "$FASTA.fai" \
    "$DB"
do
    [[ -s "$file" ]] || {
        echo "ERROR: Missing or empty file: $file" >&2
        exit 1
    }
done

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Python: $PYTHON"

"$PYTHON" -c 'import torch, gffutils, pangolin; print("Imports OK"); print("PyTorch:", torch.__version__); print("CUDA available:", torch.cuda.is_available())'

echo "Preparing combined seven-tissue SNP inputs..."
"$PYTHON" "$SOURCE/scripts/prepare_targetals_variants.py"

TEST_INPUT="$OUT/variants_test.csv"
TEST_PREFIX="$OUT/pangolin_test_${SLURM_JOB_ID}"

[[ -s "$TEST_INPUT" ]] || {
    echo "ERROR: Test input was not created." >&2
    exit 1
}

echo "Running ten-SNP Pangolin test..."
echo "Settings: unmasked; +/-1000 bp; report site changes >=0.01"

"$PANGOLIN" \
    -m False \
    -d 1000 \
    -s 0.01 \
    "$TEST_INPUT" \
    "$FASTA" \
    "$DB" \
    "$TEST_PREFIX"

[[ -s "${TEST_PREFIX}.csv" ]] || {
    echo "ERROR: Expected Pangolin output missing: ${TEST_PREFIX}.csv" >&2
    exit 1
}

echo "Test output preview:"
head -n 11 "${TEST_PREFIX}.csv"

echo "Finished: $(date)"
echo "Mapping audit: $OUT/sqtl_pairs_mapping.tsv"
echo "Full variant input: $OUT/variants.csv"
echo "Test predictions: ${TEST_PREFIX}.csv"
