#!/bin/bash
#SBATCH --job-name=pangolin_full
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=23:30:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/full_%A_%a.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/full_%A_%a.err

set -euo pipefail

: "${RUN_DIR:?Submit with RUN_DIR exported}"
: "${SLURM_ARRAY_TASK_ID:?This script requires a Slurm array}"

PYTHON="$HOME/donglab/pipelines/modules/miniconda3/bin/python"
PANGOLIN="$HOME/donglab/pipelines/modules/miniconda3/bin/pangolin"
FASTA="$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
DB="$HOME/donglab/pipelines/modules/Pangolin/references/gencode.v49.all_transcripts.db"

CHUNK_SIZE=1000
TASK=$(printf '%05d' "$SLURM_ARRAY_TASK_ID")
INPUT="$RUN_DIR/inputs/chunk_${TASK}.csv"
PREFIX="$RUN_DIR/predictions/chunk_${TASK}"
TEMP_PREFIX="$RUN_DIR/predictions/chunk_${TASK}.running_${SLURM_JOB_ID}"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"

mkdir -p "$RUN_DIR/inputs" "$RUN_DIR/predictions" "$RUN_DIR/status"

for file in "$PYTHON" "$PANGOLIN" "$FASTA" "$FASTA.fai" "$DB" "$RUN_DIR/variants.csv"; do
    [[ -s "$file" ]] || {
        echo "ERROR: Missing or empty file: $file" >&2
        exit 1
    }
done

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Task: $TASK"

"$PYTHON" - <<'PY'
import torch
print("PyTorch:", torch.__version__)
if not torch.cuda.is_available():
    raise SystemExit("ERROR: GPU unavailable; refusing CPU fallback.")
print("GPU:", torch.cuda.get_device_name(0))
PY

# Select this task's records from the frozen input.
"$PYTHON" - "$RUN_DIR/variants.csv" "$INPUT" \
    "$SLURM_ARRAY_TASK_ID" "$CHUNK_SIZE" <<'PY'
import csv
import sys
from itertools import islice

source, destination, task, size = sys.argv[1:]
start = int(task) * int(size)

with open(source, newline="") as handle:
    reader = csv.DictReader(handle)
    records = list(islice(reader, start, start + int(size)))

if not records:
    raise SystemExit("ERROR: Array task has no input records.")

with open(destination, "w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=["CHROM", "POS", "REF", "ALT"]
    )
    writer.writeheader()
    writer.writerows(records)

print(f"Selected {len(records)} variants", flush=True)
PY

echo "Scoring variants..."
"$PANGOLIN" \
    -m False \
    -d 1000 \
    -s 0.01 \
    "$INPUT" \
    "$FASTA" \
    "$DB" \
    "$TEMP_PREFIX"

# Pangolin's CSV can contain unquoted commas in its result field.
# Validate record identity using only the first four comma-separated fields.
"$PYTHON" - "$INPUT" "${TEMP_PREFIX}.csv" \
    "$RUN_DIR/status/chunk_${TASK}.summary.tsv" <<'PY'
import csv
import sys

input_file, output_file, summary_file = sys.argv[1:]
with open(input_file, newline="") as handle:
    expected = [
        (r["CHROM"], r["POS"], r["REF"], r["ALT"])
        for r in csv.DictReader(handle)
    ]

observed = []
blank = 0
with open(output_file) as handle:
    header = handle.readline().strip()
    if header != "CHROM,POS,REF,ALT,Pangolin":
        raise SystemExit(f"ERROR: Unexpected output header: {header}")
    for line in handle:
        fields = line.rstrip("\r\n").split(",", 4)
        if len(fields) != 5:
            raise SystemExit("ERROR: Malformed output record")
        observed.append(tuple(fields[:4]))
        blank += not bool(fields[4].strip())

if observed != expected:
    raise SystemExit("ERROR: Output records do not match input records.")

with open(summary_file, "w") as handle:
    handle.write("input_variants\toutput_variants\tblank_predictions\n")
    handle.write(f"{len(expected)}\t{len(observed)}\t{blank}\n")

print(f"Validated {len(observed)} records; {blank} blank predictions.")
print("Blank predictions require review of skipped-variant warnings.")
PY

mv "${TEMP_PREFIX}.csv" "${PREFIX}.csv"
date > "$RUN_DIR/status/chunk_${TASK}.complete"

echo "Output: ${PREFIX}.csv"
echo "Finished: $(date)"
