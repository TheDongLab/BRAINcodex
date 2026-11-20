#!/bin/bash
#SBATCH --job-name=circ_splicing
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/circ_exon_interval_length.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/circ_exon_interval_length.err
set -euo pipefail
module purge
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
module load deepTools
export TARGET_ROOT="${TARGET_ROOT:-$HOME/donglab/data/target_ALS}"
export ANSWER_ROOT="${ANSWER_ROOT:-$HOME/donglab/data/answer_ALS}"
export OUTDIR="${OUTDIR:-$TARGET_ROOT/QTL/circ_splicing_comparison}"
export MIN_CIRC_READS="${MIN_CIRC_READS:-2}"
export MIN_JUNCTION_READS="${MIN_JUNCTION_READS:-10}"
python - <<'PYTHON'
import csv, os, re, math, statistics
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

out=Path(os.environ['OUTDIR']); out.mkdir(parents=True,exist_ok=True)
mincirc=int(os.environ['MIN_CIRC_READS']); minj=int(os.environ['MIN_JUNCTION_READS'])
pattern=re.compile(r'(chr[^:\s]+):(\d+)-(\d+):([+-])')
roots=[('target_ALS',Path(os.environ['TARGET_ROOT'])),('answer_ALS',Path(os.environ['ANSWER_ROOT']))]
events={'unc13a_128':('chr19',17642541,17642844,'-'),'unc13a_178':('chr19',17642591,17642844,'-'),'unc13a_skip':('chr19',17641556,17642844,'-'),'stmn2_cryptic':('chr8',79611214,79616821,'-'),'stmn2_canonical':('chr8',79611214,79636801,'-')}
def write(name,rows,fields):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
issues=[]; rows=[]; stats=[]
for dataset,root in roots:
  for processed in sorted(root.glob('*/RNAseq/Processed')):
    tissue=processed.parents[1].name
    if dataset=='answer_ALS' and tissue!='iPSC_Motor_Neuron': continue
    sigfile=root/tissue/'cQTL'/'results'/f'{tissue}_cQTL.FDR0.05.txt'
    sig=set(); sig_available=sigfile.is_file()
    if sig_available:
      with sigfile.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
          hits={m.group(0) for value in row.values() if isinstance(value,str) for m in pattern.finditer(value)}
          if len(hits)==1: sig.update(hits)
          else: issues.append([str(sigfile),'row has missing or ambiguous circRNA coordinates'])
      if not sig: issues.append([str(sigfile),'no circRNA IDs parsed; cQTL subset unavailable']); sig_available=False
    else: issues.append([str(sigfile),'significant cQTL file missing'])
    for sd in sorted(x for x in processed.iterdir() if x.is_dir()):
      row={'dataset':dataset,'tissue':tissue,'sample':sd.name}; junctions=list((sd/'leafcutter'/'juncs').glob('*.junc')); counts=Counter(); jstatus='ok'
      if len(junctions)!=1: jstatus='missing' if not junctions else 'multiple_files'
      else:
        try:
          with junctions[0].open() as f:
            for line in f:
              a=line.split()
              if len(a)<6 or a[0] not in ('chr8','chr19'): continue
              try: key=(a[0],int(a[1]),int(a[2]),a[5]); n=int(a[4])
              except ValueError: continue
              for name,coord in events.items():
                if key==coord: counts[name]+=n
        except OSError as e: jstatus='read_error';issues.append([str(junctions[0]),str(e)])
      row['junction_status']=jstatus
      for name in events: row[name+'_reads']=counts[name] if jstatus=='ok' else ''
      if jstatus=='ok':
        for gene,crypt,canonical in [('unc13a',counts['unc13a_128']+counts['unc13a_178'],counts['unc13a_skip']),('stmn2',counts['stmn2_cryptic'],counts['stmn2_canonical'])]:
          depth=crypt+canonical;row[gene+'_junction_depth']=depth;row[gene+'_cryptic_fraction']=crypt/depth if depth>=minj and depth>0 else ''
      qfile=sd/'circularRNA_known_circ_percentage.txt'; circ={};qstatus='ok'
      try:
        with qfile.open() as f:
          reader=csv.DictReader(f,delimiter='\t'); required={'chrom','start','end','strand','circ_reads'}
          if not required.issubset(reader.fieldnames or []): raise ValueError('missing columns: '+','.join(sorted(required-set(reader.fieldnames or []))))
          duplicates=set()
          for item in reader:
            try:
              chrom=item['chrom'];chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
              start=int(float(item['start']));end=int(float(item['end']));strand=item['strand'].strip();reads=float(item['circ_reads'])
              if not math.isfinite(reads) or start>=end or strand not in ('+','-') or reads<0: continue
              key=f'{chrom}:{start}-{end}:{strand}'
              if key in circ: duplicates.add(key)
              circ[key]=(end-start,reads)
            except (ValueError,TypeError): continue
          for key in duplicates: circ.pop(key,None)
          if duplicates: issues.append([str(qfile),f'{len(duplicates)} duplicated IDs excluded'])
      except (OSError,ValueError) as e: qstatus='unavailable';issues.append([str(qfile),str(e)])
      row['circ_status']=qstatus
      if qstatus=='ok':
        for subset in ['all','cqtl']:
          if subset=='cqtl' and not sig_available: continue
          lengths=[v[0] for k,v in circ.items() if v[1]>=mincirc and (subset=='all' or k in sig)]
          row[subset+'_n_circs']=len(lengths);row[subset+'_median_span_bp']=statistics.median(lengths) if lengths else ''
      rows.append(row)
fields=['dataset','tissue','sample','junction_status']+[k+'_reads' for k in events]+['unc13a_junction_depth','unc13a_cryptic_fraction','stmn2_junction_depth','stmn2_cryptic_fraction','circ_status','all_n_circs','all_median_span_bp','cqtl_n_circs','cqtl_median_span_bp']
write('sample_circ_splicing.tsv',rows,fields)
for dataset,tissue in sorted({(r['dataset'],r['tissue']) for r in rows}):
  fig,axes=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
  for i,gene in enumerate(['unc13a','stmn2']):
    for j,subset in enumerate(['all','cqtl']):
      pairs=[(float(r[subset+'_median_span_bp']),float(r[gene+'_cryptic_fraction'])) for r in rows if r['dataset']==dataset and r['tissue']==tissue and r.get(subset+'_median_span_bp','')!='' and r.get(gene+'_cryptic_fraction','')!='']
      n=len(pairs);rho=pval=''
      if n>=3 and len({x for x,y in pairs})>1 and len({y for x,y in pairs})>1:
        res=spearmanr([x for x,y in pairs],[y for x,y in pairs]);rho=float(res.statistic);pval=float(res.pvalue)
      stats.append(dict(dataset=dataset,tissue=tissue,gene=gene,subset=subset,n=n,spearman_rho=rho,p_value=pval))
      ax=axes[i,j]
      if pairs: ax.scatter([x for x,y in pairs],[y for x,y in pairs],alpha=.55,s=14)
      ax.set_xlabel('Median detected circRNA span (bp)');ax.set_ylabel(gene.upper()+' cryptic junction fraction');ax.set_title(f'{subset} | n={n}');ax.set_ylim(-.03,1.03)
  fig.suptitle(dataset+' | '+tissue);fig.savefig(out/f'{dataset}_{tissue}_circ_splicing.png',dpi=160);plt.close(fig)
write('correlations.tsv',stats,['dataset','tissue','gene','subset','n','spearman_rho','p_value'])
with (out/'issues.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['file','issue']);w.writerows(issues)
print('Output:',out,'samples:',len(rows),'issues:',len(issues),'minimum circ reads:',mincirc,'minimum competing junction reads:',minj)
PYTHON
