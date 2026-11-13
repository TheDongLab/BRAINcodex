#!/bin/bash
#SBATCH --job-name=pangolin_junctions
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/junctions_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/pangolin/logs/junctions_%j.err

set -euo pipefail

"$HOME/donglab/pipelines/modules/miniconda3/bin/python" - <<'PY'
import csv
import re
from pathlib import Path
from collections import Counter
from functools import lru_cache
from contextlib import ExitStack
import gffutils

base = Path.home() / "donglab"
folder = base / "data/target_ALS/QTL/pangolin/full_20261112_163519/merged"
db = gffutils.FeatureDB(str(
    base / "pipelines/modules/Pangolin/references/"
           "gencode.v49.all_transcripts.db"
))
source = folder / "all_tissues_sqtl_pangolin.tsv"
out = folder / "junction_matched"
out.mkdir(exist_ok=True)

@lru_cache(maxsize=None)
def gene_info(gene):
    try:
        feature = db[gene]
        return feature.seqid, feature.strand
    except KeyError:
        return None

@lru_cache(maxsize=20000)
def parse_predictions(raw):
    sites = []
    for bundle in raw.split(","):
        tokens = bundle.split("|")
        gene = tokens[0]
        for token in tokens[1:]:
            match = re.fullmatch(
                r"(-?\d+):([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
                r"(?:[eE][-+]?\d+)?)", token
            )
            if match:
                sites.append((gene, int(match[1]), float(match[2])))
    return sites

def norm_chrom(chrom):
    return chrom.removeprefix("chr").upper()

extra = [
    "donor_position", "acceptor_position",
    "donor_in_scoring_window", "acceptor_in_scoring_window",
    "junction_prediction_status",
    "donor_reported_changes", "acceptor_reported_changes",
    "junction_max_abs_reported_change",
    "prediction_gene_lookup_failures",
]
counts = Counter()

with source.open(newline="") as handle, ExitStack() as stack:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames + extra

    combined = csv.DictWriter(
        stack.enter_context(
            (out / "all_tissues_junction_matched.tsv").open("w", newline="")
        ),
        fields, delimiter="\t",
    )
    combined.writeheader()
    writers = {}

    for row in reader:
        tissue = row["tissue"]
        if tissue not in writers:
            writers[tissue] = csv.DictWriter(
                stack.enter_context(
                    (out / f"{tissue}_junction_matched.tsv").open(
                        "w", newline=""
                    )
                ),
                fields, delimiter="\t",
            )
            writers[tissue].writeheader()

        row.update({field: "" for field in extra})
        chrom, strand, interval = row["junction_id"].split(":")
        start, end = map(int, interval.split("-"))
        if strand not in ("+", "-") or start >= end:
            raise RuntimeError(f"Invalid junction: {row['junction_id']}")

        donor, acceptor = (
            (start, end) if strand == "+" else (end, start)
        )
        row["donor_position"] = donor
        row["acceptor_position"] = acceptor
        status = row["pangolin_status"]

        if row["mapping_status"] != "ready":
            status = "not_submitted"
        else:
            pos = int(row["POS"])
            same_chrom = norm_chrom(chrom) == norm_chrom(row["CHROM"])
            donor_in = same_chrom and abs(donor - pos) <= 1000
            acceptor_in = same_chrom and abs(acceptor - pos) <= 1000
            row["donor_in_scoring_window"] = donor_in
            row["acceptor_in_scoring_window"] = acceptor_in

            if status == "not_scored":
                status = "not_scored"
            elif not donor_in and not acceptor_in:
                status = "both_sites_outside_window"
            else:
                hits = {"donor": [], "acceptor": []}
                values = []
                lookup_failures = set()
                for gene, offset, delta in parse_predictions(
                    row["pangolin_raw"]
                ):
                    info = gene_info(gene)
                    if info is None:
                        lookup_failures.add(gene)
                        continue
                    gene_chrom, gene_strand = info
                    if (
                        norm_chrom(gene_chrom) != norm_chrom(chrom)
                        or gene_strand != strand
                    ):
                        continue
                    site = pos + offset
                    for label, boundary in (
                        ("donor", donor), ("acceptor", acceptor)
                    ):
                        if site == boundary:
                            hits[label].append(f"{gene}:{delta}")
                            values.append(abs(delta))

                row["donor_reported_changes"] = ";".join(
                    sorted(set(hits["donor"]))
                )
                row["acceptor_reported_changes"] = ";".join(
                    sorted(set(hits["acceptor"]))
                )
                row["prediction_gene_lookup_failures"] = ";".join(
                    sorted(lookup_failures)
                )
                if values:
                    row["junction_max_abs_reported_change"] = max(values)
                    status = (
                        "both_sites_reported_change"
                        if hits["donor"] and hits["acceptor"]
                        else "donor_reported_change"
                        if hits["donor"]
                        else "acceptor_reported_change"
                    )
                else:
                    status = (
                        "both_sites_assessed_no_reported_change"
                        if donor_in and acceptor_in
                        else "one_site_assessed_no_reported_change"
                    )
                    if lookup_failures:
                        status = "gene_lookup_incomplete"

        row["junction_prediction_status"] = status
        counts[(tissue, status)] += 1
        combined.writerow(row)
        writers[tissue].writerow(row)

with (out / "junction_summary.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["tissue", "junction_prediction_status",
                     "association_rows"])
    for (tissue, status), count in sorted(counts.items()):
        writer.writerow([tissue, status, count])
        print(f"{tissue}\t{status}\t{count:,}")

print(f"\nOutputs: {out}")
print("Reported changes use the original >=0.01 screening threshold.")
print("QTL betas are retained unchanged; no direction-concordance claim.")
PY
