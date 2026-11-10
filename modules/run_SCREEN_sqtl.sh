#!/bin/bash
#SBATCH --job-name=SCREEN_sqtl
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/SCREEN_sQTL_PSI_ranking/run_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/SCREEN_sQTL_PSI_ranking/run_%j.err
# Self-contained SCREEN overlap ranking and genotype/PSI plotting.
# Submit: sbatch ~/donglab/pipelines/scripts/QTL/run_SCREEN_sqtl.sh
# Results in ~/donglab/data/target_ALS/QTL/SCREEN_sQTL_PSI_ranking/
set -euo pipefail
TASK_PYTHON=/home/zw529/donglab/pipelines/modules/miniconda3/bin/python
mkdir -p /home/zw529/donglab/data/target_ALS/QTL/SCREEN_sQTL_PSI_ranking
"$TASK_PYTHON" -I -u - <<'SCREEN_RANK_PYTHON'
#!/usr/bin/env python3
"""Rank SCREEN-overlapping sQTLs by observed homozygote PSI difference. Coordinates: SNP 1-based; SCREEN BED 0-based half-open."""
import argparse
import csv
from bisect import bisect_left
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

TISSUES = [
    "Cerebellum", "Cervical_Spinal_Cord", "Cortex", "Frontal_Cortex",
    "Lumbar_Spinal_Cord", "Motor_Cortex", "Spinal_Cord",
]


def matrix(path, wanted):
    data = {}
    with path.open() as handle:
        samples = handle.readline().split()[1:]
        if len(samples) != len(set(samples)):
            raise ValueError(f"Duplicate subjects: {path}")
        for line in handle:
            key, _, tail = line.rstrip("\n").partition("\t")
            if key not in wanted:
                continue
            if key in data:
                raise ValueError(f"Duplicate matrix row: {key}")
            data[key] = np.array([
                np.nan if x in ("NA", "NaN", "nan", ".", "") else float(x)
                for x in tail.split("\t")
            ])
            if len(data[key]) != len(samples):
                raise ValueError(f"Wrong row length: {key}")
    return samples, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path.home() / "donglab/data/target_ALS")
    parser.add_argument('--screen', type=Path, default=Path('/home/zw529/donglab/references/genome/Homo_sapiens/ENCODE/SCREEN'))
    args = parser.parse_args()
    root = args.root
    out = root / "QTL/SCREEN_sQTL_PSI_ranking"
    screen = args.screen
    overlap_counts = {}
    missing_locations = []
    sigs, locations, targets = {}, {}, set()
    for tissue in TISSUES:
        q = root / tissue / "sQTL"
        s = pd.read_csv(q / "results" / f"{tissue}_sQTL.FDR0.05.txt",
                        sep="\t", dtype={"junction_id": str, "snpid": str})
        s = s[pd.to_numeric(s.FDR, errors="raise") <= 0.05]
        s = s.drop_duplicates(["junction_id", "snpid"])
        sigs[tissue] = s
        wanted, loc = set(s.snpid), {}
        with (q / "snp_location.txt").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["snpid"] not in wanted:
                    continue
                coord = (row["chr"].removeprefix("chr"), int(row["pos"]))
                if row["snpid"] in loc and loc[row["snpid"]] != coord:
                    raise ValueError(f"Ambiguous SNP location: {tissue} {row['snpid']}")
                loc[row["snpid"]] = coord
                targets.add(coord)
        locations[tissue] = loc

    # Index only candidate positions; stream the complete registry once.
    by_chr = defaultdict(list)
    for chrom, pos in targets:
        if pos < 1:
            raise ValueError('SNP positions must be 1-based positive integers')
        by_chr[chrom].append(pos - 1)
    by_chr = {chrom: sorted(set(pos)) for chrom, pos in by_chr.items()}
    overlaps = defaultdict(set)
    for cls in ('PLS', 'pELS', 'dELS'):
        with (screen / f'GRCh38-cCREs.{cls}.bed').open() as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.rstrip().split('\t')
                chrom, lo, hi = fields[0].removeprefix('chr'), int(fields[1]), int(fields[2])
                if lo < 0 or hi <= lo or fields[-1] != cls:
                    raise ValueError(f'Invalid SCREEN record: {line}')
                # v4 has both a region ID (col 4) and an EH38E cCRE ID (col 5).
                element = fields[4]
                positions = by_chr.get(chrom, [])
                for pos0 in positions[bisect_left(positions, lo):bisect_left(positions, hi)]:
                    overlaps[chrom, pos0 + 1].add((cls, element, f'chr{chrom}:{lo}-{hi}'))
    out.mkdir(parents=True, exist_ok=True)
    overlap_rows = []
    for tissue in TISSUES:
        s = sigs[tissue]
        keep = []
        for row in s.to_dict('records'):
            coord = locations[tissue].get(row['snpid'])
            hits = overlaps.get(coord, set())
            keep.append(bool(hits))
            if coord is None:
                missing_locations.append(dict(row, tissue=tissue, reason='missing_snp_location'))
            for cls, element, interval in sorted(hits):
                overlap_rows.append(dict(tissue=tissue, junction_id=row['junction_id'], snpid=row['snpid'],
                                         regulatory_class=cls, cCRE_id=element, interval_bed=interval))
        sigs[tissue] = s.loc[keep].copy()
        overlap_counts[tissue] = (len(s), len(sigs[tissue]))
    pd.DataFrame(overlap_rows, columns=['tissue','junction_id','snpid','regulatory_class','cCRE_id','interval_bed']).to_csv(out/'SCREEN_overlaps.tsv', sep='\t', index=False)
    pd.DataFrame(missing_locations, columns=['tissue','junction_id','snpid','reason']).to_csv(out/'missing_locations.tsv', sep='\t', index=False)

    alleles = defaultdict(set)
    with (root / "QTL/plink/joint_all_chrs_matrixEQTL.raw").open() as handle:
        for token in handle.readline().split()[6:]:
            m = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([^:]+):([^_]+)_(.+)", token)
            if m:
                chrom, pos, ref, alt, counted = m.groups()
                coord = (chrom, int(pos))
                if coord in targets:
                    alleles[coord].add((ref, alt, counted))

    out.mkdir(parents=True, exist_ok=True)
    all_top, audit = [], []
    metrics = [
        "tissue", "regulatory_classes", "cCRE_ids", "cCRE_intervals_bed", "variant_chr", "variant_pos", "ref", "alt",
        "plink_counted_allele", "n_ref_ref", "n_het", "n_hom_alt",
        "mean_ref_ref_pct", "mean_het_pct", "mean_hom_alt_pct",
        "delta_alt_minus_ref_pp", "abs_delta_pp",
    ]
    display = [
        "top_rank", "regulatory_classes", "cCRE_ids", "cCRE_intervals_bed", "junction_id", "snpid", "ref", "alt",
        "n_ref_ref", "n_hom_alt", "mean_ref_ref_pct", "mean_hom_alt_pct",
        "delta_alt_minus_ref_pp", "abs_delta_pp", "FDR",
    ]
    for tissue in TISSUES:
        print(f"\nLoading {tissue}...", flush=True)
        q, s = root / tissue / "sQTL", sigs[tissue]
        records, excluded = [], []
        if len(s):
            subjects, psi = matrix(q / f"splicing_{tissue}.txt", set(s.junction_id))
            gs, gts = matrix(q / f"snp_{tissue}.txt", set(s.snpid))
            if set(subjects) != set(gs):
                raise ValueError(f"PSI/genotype subject sets differ: {tissue}")
            gi = {sample: i for i, sample in enumerate(gs)}
            order = [gi[sample] for sample in subjects]
            for row in s.to_dict("records"):
                j, snp = row["junction_id"], row["snpid"]
                coord = locations[tissue].get(snp)
                aa = alleles.get(coord, set())
                reason = None
                if j not in psi or snp not in gts:
                    reason = "missing_matrix_row"
                elif len(aa) != 1:
                    reason = "missing_or_ambiguous_alleles"
                else:
                    ref, alt, counted = next(iter(aa))
                    if counted not in (ref, alt) or "," in alt or ref == alt:
                        reason = "unsupported_alleles"
                if reason:
                    excluded.append(dict(row, tissue=tissue, reason=reason))
                    continue
                x, g = psi[j], gts[snp][order]
                if np.any(np.isfinite(x) & ((x < 0) | (x > 1))):
                    raise ValueError(f"Non-PSI value: {tissue} {j}")
                if np.any(np.isfinite(g) & ~np.isin(g, [0, 1, 2])):
                    raise ValueError(f"Non-hardcall genotype: {tissue} {snp}")
                altg = 2 - g if counted == ref else g
                hits = sorted(overlaps[coord])
                r = dict(row, tissue=tissue,
                         regulatory_classes=';'.join(sorted({h[0] for h in hits})),
                         cCRE_ids=';'.join(h[1] for h in hits),
                         cCRE_intervals_bed=';'.join(h[2] for h in hits),
                         variant_chr="chr" + coord[0], variant_pos=coord[1],
                         ref=ref, alt=alt, plink_counted_allele=counted)
                for value, label in [(0, "ref_ref"), (1, "het"), (2, "hom_alt")]:
                    v = x[(altg == value) & np.isfinite(x)]
                    r["n_" + label] = len(v)
                    r["mean_" + label + "_pct"] = 100 * float(v.mean()) if len(v) else np.nan
                r["delta_alt_minus_ref_pp"] = r["mean_hom_alt_pct"] - r["mean_ref_ref_pct"]
                r["abs_delta_pp"] = abs(r["delta_alt_minus_ref_pp"])
                if r["n_ref_ref"] < 5 or r["n_hom_alt"] < 5:
                    excluded.append(dict(r, reason="fewer_than_5_usable_subjects_in_one_or_both_homozygote_groups"))
                else:
                    records.append(r)
        columns = list(dict.fromkeys(list(s.columns) + metrics))
        ranked = pd.DataFrame(records, columns=columns)
        ranked = ranked.sort_values(
            ["abs_delta_pp", "FDR", "junction_id", "snpid"],
            ascending=[False, True, True, True],
        ).reset_index(drop=True)
        ranked.insert(0, "rank", np.arange(1, len(ranked) + 1))
        top = ranked.loc[ranked.groupby("junction_id", sort=False).cumcount() < 2].head(100).copy()
        top.insert(0, "top_rank", np.arange(1, len(top) + 1))
        ranked.to_csv(out / f"{tissue}.ranked_qualifying.tsv", sep="\t", index=False, na_rep="NA")
        top.to_csv(out / f"{tissue}.top100.tsv", sep="\t", index=False, na_rep="NA")
        pd.DataFrame(excluded, columns=columns + ["reason"]).to_csv(
            out / f"{tissue}.excluded.tsv", sep="\t", index=False, na_rep="NA")
        print(f"\n{tissue}: {len(s):,} SCREEN-overlapping pairs; {len(ranked):,} qualifying; {len(top)} selected")
        print(top[display].to_string(index=False, float_format=lambda v: f"{v:.6g}")
              if len(top) else "No qualifying pairs.", flush=True)
        plot_rows = []
        for r in top.to_dict('records'):
            x, g = psi[r['junction_id']], gts[r['snpid']][order]
            altg = 2-g if r['plink_counted_allele']==r['ref'] else g
            for sample, value, dosage in zip(subjects, x, altg):
                if np.isfinite(value) and np.isfinite(dosage):
                    plot_rows.append(dict(top_rank=r['top_rank'], junction_id=r['junction_id'], snpid=r['snpid'],
                                          subject=sample, PSI_pct=100*value, ALT_dosage=int(dosage)))
        pd.DataFrame(plot_rows, columns=['top_rank','junction_id','snpid','subject','PSI_pct','ALT_dosage']).to_csv(out/f'{tissue}.top100_subjects.tsv', sep='\t', index=False)
        all_top.append(top)
        audit.append(dict(tissue=tissue, significant_pairs=overlap_counts[tissue][0], SCREEN_overlapping_pairs=len(s), no_SCREEN_overlap_pairs=overlap_counts[tissue][0]-len(s)-sum(r["tissue"]==tissue for r in missing_locations), missing_location_pairs=sum(r["tissue"]==tissue for r in missing_locations),
                          qualifying_pairs=len(ranked), selected_pairs=len(top),
                          excluded_pairs=len(excluded)))
        if len(s):
            del psi, gts
    pd.concat(all_top, ignore_index=True).to_csv(
        out / "top100_per_tissue.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(audit).to_csv(out / "audit.tsv", sep="\t", index=False)
    print(f"\nSaved results to {out}", flush=True)


if __name__ == "__main__":
    main()

SCREEN_RANK_PYTHON
"$TASK_PYTHON" -I -u - <<'SCREEN_PLOT_PYTHON'
#!/usr/bin/env python3
"""Plot observed PSI by REF/ALT genotype for the SCREEN-selected top pairs."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir', type=Path, default=Path('/home/zw529/donglab/data/target_ALS/QTL/SCREEN_sQTL_PSI_ranking'))
    a = ap.parse_args()
    dest = a.outdir / 'plots'
    dest.mkdir(exist_ok=True)
    for path in sorted(a.outdir.glob('*.top100.tsv')):
        tissue = path.name.removesuffix('.top100.tsv')
        top = pd.read_csv(path, sep='\t')
        subjects = pd.read_csv(a.outdir / f'{tissue}.top100_subjects.tsv', sep='\t')
        if top.empty:
            continue
        nrows = (len(top)+1)//2
        fig, axes = plt.subplots(nrows, 2, figsize=(13, 3.8*nrows), squeeze=False)
        rng = np.random.default_rng(42)
        for ax, row in zip(axes.flat, top.to_dict('records')):
            sub = subjects[subjects.top_rank == row['top_rank']]
            groups = [sub.loc[sub.ALT_dosage == g, 'PSI_pct'].to_numpy() for g in (0,1,2)]
            for g, (values, color) in enumerate(zip(groups, ['#3979ad','#9a9a9a','#cb6743'])):
                if len(values):
                    ax.scatter(g+rng.uniform(-.12,.12,len(values)), values, s=11, alpha=.45, color=color, linewidths=0)
                    ax.hlines(values.mean(), g-.23, g+.23, color='black', linewidth=2)
            ax.set_xticks([0,1,2])
            ax.set_xticklabels([f'{row["ref"]}/{row["ref"]}\nn={len(groups[0])}',f'{row["ref"]}/{row["alt"]}\nn={len(groups[1])}',f'{row["alt"]}/{row["alt"]}\nn={len(groups[2])}'])
            ax.set_ylim(-3,103)
            ax.set_ylabel('Observed PSI (%)')
            ax.set_title(f'#{row["top_rank"]}  {row["snpid"]}  |  {row["regulatory_classes"]}\n{row["junction_id"]}\nALT−REF = {row["delta_alt_minus_ref_pp"]:+.1f} pp; QTL FDR = {row["FDR"]:.2g}',fontsize=9)
            ax.spines[['top','right']].set_visible(False)
        for ax in list(axes.flat)[len(top):]:
            ax.set_visible(False)
        fig.suptitle(f'{tissue}: SCREEN-overlapping sQTLs\nObserved PSI; black lines show means. No brain-specific SCREEN filtering.', fontsize=12)
        fig.tight_layout(rect=(0,0,1,.96))
        for ext in ('png','pdf'):
            fig.savefig(dest/f'{tissue}.top100.{ext}',dpi=180,bbox_inches='tight')
        plt.close(fig)
        print(f'Plotted {tissue}: {len(top)} pairs', flush=True)

if __name__ == '__main__':
    main()

SCREEN_PLOT_PYTHON
