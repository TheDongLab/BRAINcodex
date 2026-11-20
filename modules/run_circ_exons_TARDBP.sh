#!/bin/bash
#SBATCH --job-name=circ_exons_TARDBP
#SBATCH --time=08:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/circ_exons_TARDBP.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/circ_exons_TARDBP.err

set -euo pipefail
# Submit with: sbatch run_circ_exons_TARDBP.sh
# Optional overrides: export TARGET_ROOT=/path/to/target_ALS GTF=/path/to/annotation.gtf OUTDIR=/path/to/output
TARGET_ROOT="${TARGET_ROOT:-$HOME/donglab/data/target_ALS}"
ANSWER_ROOT="${ANSWER_ROOT:-$HOME/donglab/data/answer_ALS}"
GTF="${GTF:-$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf}"
OUTDIR="${OUTDIR:-$TARGET_ROOT/QTL/circ_exon_TARDBP_summary}"
python - --root "$TARGET_ROOT" --answer-root "$ANSWER_ROOT" --gtf "$GTF" --out "$OUTDIR" <<'PYTHON_ANALYSIS'
#!/usr/bin/env python3
"""Summarize significant circRNA QTL intervals and TARDBP sample inclusion."""
import argparse
import csv
import gzip
import re
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, default=Path.home() / 'donglab/data/target_ALS')
p.add_argument('--answer-root', type=Path, default=Path.home() / 'donglab/data/answer_ALS')
p.add_argument('--gtf', type=Path, default=Path.home() / 'donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf')
p.add_argument('--out', type=Path, default=Path.home() / 'donglab/data/target_ALS/QTL/circ_exon_TARDBP_summary')
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

def write(name, fields, rows):
    with (a.out / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

event_re = re.compile(r'(chr[^:\s]+):(\d+)-(\d+):([+-])')
events = []
issues = []
for dataset, root in [('target_ALS', a.root), ('answer_ALS', a.answer_root / 'iPSC_Motor_Neuron')]:
  for kind, label in [('cQTL', 'cQTL'), ('circ_eQTL', 'eQTL')]:
    pattern = f'*/{kind}/results/*_{label}.FDR0.05.txt' if dataset == 'target_ALS' else f'{kind}/results/*_{label}.FDR0.05.txt'
    for file in sorted(root.glob(pattern)):
        tissue = file.parents[2].name
        with file.open() as f:
            reader = csv.DictReader(f, delimiter='\t')
            seen = set()
            for row in reader:
                candidates = [(k, event_re.search(str(v or ''))) for k, v in row.items() if k]
                matches = [(k, m) for k, m in candidates if m]
                if not matches:
                    issues.append([str(file), 'row lacks chr:start-end:strand event']); continue
                coords = {m.group(0) for _, m in matches}
                if len(coords) != 1:
                    issues.append([str(file), 'row has ambiguous event coordinates']); continue
                event = next(iter(coords))
                if event in seen: continue
                seen.add(event)
                m = event_re.fullmatch(event)
                events.append(dict(dataset=dataset, tissue=tissue, qtl_type=kind, event=event, chrom=m[1], start=int(m[2]), end=int(m[3]), strand=m[4], source=str(file)))

# Count distinct annotated exon intervals fully within the genomic circRNA span,
# on its strand. Coordinates are treated as 1-based inclusive, like GTF.
exons = defaultdict(set)
opener = gzip.open if str(a.gtf).endswith('.gz') else open
with opener(a.gtf, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        v = line.rstrip('\n').split('\t')
        if len(v) >= 7 and v[2] == 'exon':
            exons[(v[0], v[6])].add((int(v[3]), int(v[4])))
for key in exons: exons[key] = sorted(exons[key])
for e in events:
    e['exon_count'] = sum(1 for s, t in exons[(e['chrom'], e['strand'])] if s >= e['start'] and t <= e['end'])
write('circ_exon_counts.tsv', ['dataset','tissue','qtl_type','event','chrom','start','end','strand','exon_count','source'], events)
for dataset, tissue in sorted({(e['dataset'], e['tissue']) for e in events}):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, kind in zip(axes, ['cQTL','circ_eQTL']):
        vals = [e['exon_count'] for e in events if e['dataset'] == dataset and e['tissue'] == tissue and e['qtl_type'] == kind]
        counts = Counter(vals)
        if counts: ax.bar(sorted(counts), [counts[k] for k in sorted(counts)], width=.85)
        else: ax.text(.5, .5, 'No events found', ha='center', transform=ax.transAxes)
        ax.set_title(f'{kind} (n={len(vals)} circRNAs)'); ax.set_xlabel('Distinct annotated exons fully within circRNA span'); ax.set_ylabel('Number of circRNAs')
    fig.suptitle(f'{dataset} | {tissue}'); fig.savefig(a.out / f'{dataset}_{tissue}_circ_exon_histogram.png', dpi=180); plt.close(fig)

samples = []
for dataset, root in [('target_ALS', a.root), ('answer_ALS', a.answer_root)]:
  processed_dirs = sorted(root.glob('*/RNAseq/Processed'))
  if dataset == 'answer_ALS': processed_dirs = [x for x in processed_dirs if x.parents[1].name == 'iPSC_Motor_Neuron']
  for processed in processed_dirs:
    tissue = processed.parents[1].name
    matrix = root / tissue / 'eQTL' / f'expression_{tissue}.txt'
    if not matrix.exists():
        issues.append([str(matrix), 'eQTL expression matrix missing']); continue
    with matrix.open() as f: header = f.readline().rstrip('\r\n').split('\t')
    if len(header) < 2:
        issues.append([str(matrix), 'matrix header is not tab-delimited or has fewer than two columns']); continue
    eqtl = set(header[1:]); dirs = {x.name for x in processed.iterdir() if x.is_dir()}
    overlap = dirs & eqtl
    if not overlap:
        issues.append([str(matrix), 'no exact matches between matrix header and sample directories; membership not assigned']); continue
    for sample in sorted(dirs):
        path = processed / sample / 'normalization.tab'
        row = dict(dataset=dataset, tissue=tissue, sample=sample, in_eqtl='yes' if sample in eqtl else 'no', status='', TPM='', FPKM='', raw_count='')
        try:
            with path.open() as f:
                found = [r for r in csv.DictReader(f, delimiter='\t') if r.get('gene_id','').split('.')[0] == 'ENSG00000120948']
            if len(found) != 1: row['status'] = f'TARDBP_rows_{len(found)}'
            else:
                row.update({k: found[0].get(k, '') for k in ('TPM','FPKM','raw_count')}); row['status'] = 'ok'
        except OSError as exc: row['status'] = 'file_error'; issues.append([str(path), str(exc)])
        samples.append(row)
    for sample in sorted(eqtl - dirs): issues.append([str(matrix), f'eQTL sample has no matching processed directory: {sample}'])
write('TARDBP_sample_expression.tsv', ['dataset','tissue','sample','in_eqtl','status','TPM','FPKM','raw_count'], samples)
summary = []
for dataset, tissue in sorted({(r['dataset'], r['tissue']) for r in samples}):
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    vals = []
    for group in ['yes','no']:
        data = []
        for r in samples:
            if r['dataset'] == dataset and r['tissue'] == tissue and r['in_eqtl'] == group and r['status'] == 'ok':
                try: data.append(float(r['TPM']))
                except ValueError: pass
        data.sort(); vals.append(data)
        median = '' if not data else (data[(len(data)-1)//2] + data[len(data)//2])/2
        summary.append(dict(dataset=dataset, tissue=tissue, in_eqtl=group, n=len(data), median_TPM=median, zero_TPM=sum(x==0 for x in data)))
    nonempty = [(i+1,v) for i,v in enumerate(vals) if v]
    if nonempty:
        ax.boxplot([v for _,v in nonempty], positions=[i for i,_ in nonempty], showfliers=True)
    ax.set_xticks([1,2], ['Included in eQTL','Other RNA-seq']); ax.set_ylabel('TARDBP TPM'); ax.set_title(f'{dataset} | {tissue}')
    fig.savefig(a.out / f'{dataset}_{tissue}_TARDBP_TPM.png', dpi=180); plt.close(fig)
write('TARDBP_summary.tsv', ['dataset','tissue','in_eqtl','n','median_TPM','zero_TPM'], summary)
with (a.out / 'issues.tsv').open('w') as f:
    w=csv.writer(f, delimiter='\t'); w.writerow(['source','issue']); w.writerows(issues)
print(f'Output: {a.out}; {len(events)} unique tissue/QTL circRNAs; {len(samples)} RNA-seq samples; {len(issues)} issues')

PYTHON_ANALYSIS
