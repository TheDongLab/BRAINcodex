#!/bin/bash
#SBATCH --job-name=pangolin_merge
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/merge_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/merge_%j.err

set -euo pipefail

PYTHON="$HOME/donglab/pipelines/modules/miniconda3/bin/python"
RUN="$HOME/donglab/data/target_ALS/QTL/pangolin/full_20261112_163519"

"$PYTHON" - "$RUN" <<'PY'
import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path
from contextlib import ExitStack

run = Path(sys.argv[1])
out = run / "merged"
size = 1000

with (run / "variants.csv").open(newline="") as handle:
    expected = [
        (r["CHROM"], r["POS"], r["REF"], r["ALT"])
        for r in csv.DictReader(handle)
    ]

if not expected or len(set(expected)) != len(expected):
    raise SystemExit("ERROR: Input is empty or contains duplicate variants.")

nchunks = math.ceil(len(expected) / size)
missing = []
for i in range(nchunks):
    task = f"{i:05d}"
    for path in (
        run / "status" / f"chunk_{task}.complete",
        run / "predictions" / f"chunk_{task}.csv",
    ):
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(str(path))

if missing:
    raise SystemExit("ERROR: Missing/empty required files:\n" +
                     "\n".join(missing))

# Validate every chunk before generating merged outputs.
raw_predictions = {}
for i in range(nchunks):
    path = run / "predictions" / f"chunk_{i:05d}.csv"
    observed = []
    with path.open() as handle:
        if handle.readline().strip() != "CHROM,POS,REF,ALT,Pangolin":
            raise SystemExit(f"ERROR: Unexpected header: {path}")
        for line in handle:
            fields = line.rstrip("\r\n").split(",", 4)
            if len(fields) != 5:
                raise SystemExit(f"ERROR: Malformed record: {path}")
            key = tuple(fields[:4])
            observed.append(key)
            if key in raw_predictions:
                raise SystemExit(f"ERROR: Duplicate output variant: {key}")
            raw_predictions[key] = fields[4]

    if observed != expected[i * size:(i + 1) * size]:
        raise SystemExit(f"ERROR: Input/output mismatch: {path}")

print(f"Verified {nchunks} chunks and {len(expected):,} variants.",
      flush=True)

out.mkdir(exist_ok=True)
annotations = {}
counts = Counter()

fields = [
    "variant_id", "CHROM", "POS", "REF", "ALT",
    "pangolin_status", "max_abs_reported_change",
    "strongest_site_changes", "pangolin_raw",
]

with (out / "variant_predictions.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fields, delimiter="\t")
    writer.writeheader()

    for key in expected:
        chrom, pos, ref, alt = key
        raw = raw_predictions[key]
        scores = []

        for bundle in raw.split(","):
            parts = bundle.split("|")
            gene = parts[0]
            for token in parts[1:]:
                match = re.fullmatch(
                    r"(-?\d+):([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
                    r"(?:[eE][-+]?\d+)?)", token
                )
                if match:
                    offset, value = int(match[1]), float(match[2])
                    scores.append((gene, offset, value))

        status = (
            "not_scored" if not raw.strip()
            else "scored_reported_change" if scores
            else "scored_no_reported_change"
        )
        maximum = max((abs(s[2]) for s in scores), default=None)
        strongest = [
            f"{gene}|{chrom}:{int(pos)+offset}|offset={offset}|delta={value}"
            for gene, offset, value in scores
            if abs(value) == maximum
        ]

        record = dict(
            variant_id=":".join(key),
            CHROM=chrom, POS=pos, REF=ref, ALT=alt,
            pangolin_status=status,
            max_abs_reported_change=(
                maximum if maximum is not None else ""
            ),
            strongest_site_changes=";".join(strongest),
            pangolin_raw=raw,
        )
        annotations[record["variant_id"]] = record
        counts[status] += 1
        writer.writerow(record)

audit = run / "sqtl_pairs_mapping.tsv"
extra = [
    "pangolin_status", "max_abs_reported_change",
    "strongest_site_changes", "pangolin_raw",
]
tissue_counts = Counter()

with audit.open(newline="") as handle, ExitStack() as stack:
    reader = csv.DictReader(handle, delimiter="\t")
    output_fields = reader.fieldnames + extra
    combined = csv.DictWriter(
        stack.enter_context(
            (out / "all_tissues_sqtl_pangolin.tsv").open("w", newline="")
        ),
        output_fields, delimiter="\t",
    )
    combined.writeheader()
    writers = {}

    for row in reader:
        tissue = row["tissue"]
        if tissue not in writers:
            writers[tissue] = csv.DictWriter(
                stack.enter_context(
                    (out / f"{tissue}_sqtl_pangolin.tsv").open(
                        "w", newline=""
                    )
                ),
                output_fields, delimiter="\t",
            )
            writers[tissue].writeheader()

        if row["mapping_status"] == "ready":
            annotation = annotations.get(row["variant_id"])
            if annotation is None:
                raise SystemExit(
                    f"ERROR: Ready variant missing: {row['variant_id']}"
                )
            row.update({field: annotation[field] for field in extra})
        else:
            row.update({field: "" for field in extra})
            row["pangolin_status"] = (
                "not_submitted:" + row["mapping_status"]
            )

        combined.writerow(row)
        writers[tissue].writerow(row)
        tissue_counts[(tissue, row["pangolin_status"])] += 1

with (out / "tissue_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["tissue", "status", "association_rows"])
    for (tissue, status), count in sorted(tissue_counts.items()):
        writer.writerow([tissue, status, count])

print("\nUnique-variant summary:")
for status, count in sorted(counts.items()):
    print(f"  {status}: {count:,}")
print("\nAssociation-row summary:")
for (tissue, status), count in sorted(tissue_counts.items()):
    print(f"  {tissue}\t{status}\t{count:,}")
print(f"\nOutputs: {out}")
print("Original QTL betas retained unchanged.")
print("Scores are across modeled tissues, not brain-specific.")
print("Junction matching and coordinate-convention checks remain pending.")
PY
