#!/bin/bash
#SBATCH --job-name=circ_exon_bp_intervals
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
export GTF="${GTF:-$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf}"
export MIN_CIRC_READS="${MIN_CIRC_READS:-2}"
export MIN_JUNCTION_READS="${MIN_JUNCTION_READS:-10}"
python - <<'PYTHON'
import csv, os, re, math, statistics, gzip, bisect
from pathlib import Path
from collections import Counter, defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

out=Path(os.environ['OUTDIR']); out.mkdir(parents=True,exist_ok=True)
mincirc=int(os.environ['MIN_CIRC_READS']); minj=int(os.environ['MIN_JUNCTION_READS'])
pattern=re.compile(r'(chr[^:\s]+):(\d+)-(\d+):([+-])')
def normalize_id(value):
  value=str(value or '').strip()
  chrom=r'(?:chr)?([0-9]+|X|Y|M|MT)'
  patterns=[(rf'^{chrom}:(\d+)[-:](\d+):([+-])$',False), (rf'^{chrom}:([+-]):(\d+)-(\d+)$',True), (rf'^{chrom}:(\d+)[-:](\d+)$',False)]
  for pat,first in patterns:
    m=re.fullmatch(pat,value,re.I)
    if m:
      parts=m.groups(); c=parts[0].upper(); c='M' if c=='MT' else c
      if first: strand,start,end=parts[1:]
      else: start,end=parts[1:3]; strand=parts[3] if len(parts)>3 else ''
      if int(start)>=int(end): return None
      return f'chr{c}:{int(start)}-{int(end)}'+(':'+strand if strand else '')
  return None

roots=[('target_ALS',Path(os.environ['TARGET_ROOT'])),('answer_ALS',Path(os.environ['ANSWER_ROOT']))]
events={'unc13a_128':('chr19',17642541,17642844,'-'),'unc13a_178':('chr19',17642591,17642844,'-'),'unc13a_skip':('chr19',17641556,17642844,'-'),'stmn2_cryptic':('chr8',79611214,79616821,'+'),'stmn2_canonical':('chr8',79611214,79636801,'+')}
def write(name,rows,fields):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
# Exon intervals are distinct genomic features across transcripts, not a reconstructed isoform.
# Circ coordinates are assumed 1-based inclusive, matching GTF; verify the source convention.
gtf=Path(os.environ['GTF']); exon_sets=defaultdict(set)
if not gtf.is_file(): raise SystemExit('Missing GTF: '+str(gtf))
with (gzip.open(gtf,'rt') if str(gtf).endswith('.gz') else gtf.open()) as f:
  for line in f:
    if line.startswith('#'): continue
    v=line.rstrip('\n').split('\t')
    if len(v)>=7 and v[2]=='exon': exon_sets[(v[0],v[6])].add((int(v[3]),int(v[4])))
exons={k:sorted(v) for k,v in exon_sets.items()}; exon_starts={k:[a for a,b in v] for k,v in exons.items()}; exon_cache={}
def count_exons(event):
  if event not in exon_cache:
    m=pattern.fullmatch(event); chrom,start,end,strand=m.group(1),int(m.group(2)),int(m.group(3)),m.group(4)
    key=(chrom,strand); intervals=exons.get(key,[]); starts=exon_starts.get(key,[])
    left=bisect.bisect_left(starts,start); right=bisect.bisect_right(starts,end)
    exon_cache[event]=sum(b<=end for a,b in intervals[left:right])
  return exon_cache[event]
remap={}
for x in ['Motor Cortex Lateral','Motor Cortex Medial','Lateral Motor Cortex','Medial Motor Cortex','Primary Motor Cortex L','Primary Motor Cortex M','Cortex_Motor_Unspecified','Cortex_Motor_BA4','BA4 Motor Cortex','Lateral_motor_cortex','Motor_Cortex_Lateral','Motor_Cortex_Medial','Lateral_Motor_Cortex','Medial_Motor_Cortex','Primary_Motor_Cortex_L','Primary_Motor_Cortex_M','Motor_Cortex']: remap[x]='Motor_Cortex'
for x in ['Spinal_Cord_Cervical','Cervical Spinal Cord','Cervical_spinal_cord','Spinal_cord_Cervical','Cervical_Spinal_Cord']: remap[x]='Cervical_Spinal_Cord'
for x in ['Lumbar Spinal Cord','Spinal_Cord_Lumbosacral','Lumbosacral_Spinal_Cord','Lumbar_spinal_cord','Lumbar_Spinal_Cord']: remap[x]='Lumbar_Spinal_Cord'
remap.update({'Frontal Cortex':'Frontal_Cortex','Frontal_Cortex':'Frontal_Cortex','Cerebellum':'Cerebellum'})
issues=[]; rows=[]; stats=[]
for dataset,root in roots:
  for processed in sorted(root.glob('*/RNAseq/Processed')):
    original_tissue=processed.parents[1].name
    if dataset=='answer_ALS':
      if original_tissue!='iPSC_Motor_Neuron': continue
      tissue=original_tissue
    else:
      if original_tissue not in remap: continue
      tissue=remap[original_tissue]
    sigsets={}; sig_available={}
    for subset,folder,label in [('cqtl','cQTL','cQTL'),('circ_eqtl','circ_eQTL','eQTL')]:
      sigfile=root/tissue/folder/'results'/f'{tissue}_{label}.FDR0.05.txt'
      if not sigfile.is_file() and original_tissue!=tissue:
        sigfile=root/original_tissue/folder/'results'/f'{original_tissue}_{label}.FDR0.05.txt'
      sig=set(); sig_available[subset]=sigfile.is_file()
      if sig_available[subset]:
        with sigfile.open() as f:
          reader=csv.DictReader(f,delimiter='\t')
          aliases=['geneid','gene_id','gene','circ_id','circid','event_id','phenotype_id','circRNA','circRNA_ID','circ']
          cols={str(k).strip().lower():k for k in (reader.fieldnames or [])}
          col=next((cols[x.lower()] for x in aliases if x.lower() in cols),None)
          if not col:
            issues.append([str(sigfile),'no recognized circRNA ID column; columns: '+', '.join(reader.fieldnames or [])]); sig_available[subset]=False
          else:
            bad=0; examples=[]
            for record in reader:
              raw=record.get(col,''); event=normalize_id(raw)
              if event: sig.add(event)
              else:
                bad+=1
                if len(examples)<3: examples.append(str(raw))
            if bad: issues.append([str(sigfile),f'{bad} unparsed IDs in column {col}; examples: '+repr(examples)])
            if not sig: issues.append([str(sigfile),'no circRNA IDs parsed; '+subset+' subset unavailable']);sig_available[subset]=False
      else: issues.append([str(sigfile),subset+' significant results file missing'])
      sigsets[subset]=sig
    for sd in sorted(x for x in processed.iterdir() if x.is_dir()):
      row={'dataset':dataset,'tissue':tissue,'original_tissue':original_tissue,'sample':sd.name}; junctions=list((sd/'leafcutter'/'juncs').glob('*.junc')); counts=Counter(); jstatus='ok'
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
              chrom=item['chrom'].strip();chrom='chr'+re.sub(r'^chr','',chrom,flags=re.I).upper();chrom='chrM' if chrom=='chrMT' else chrom
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
        by_unstranded=defaultdict(list)
        for key in circ: by_unstranded[key.rsplit(':',1)[0]].append(key)
        for subset in ['all','cqtl','circ_eqtl']:
          if subset!='all' and not sig_available[subset]: continue
          sig=sigsets.get(subset,set())
          # Unstranded QTL IDs are used only when exactly one circRNA matches in this sample.
          def selected(k):
            unstranded=k.rsplit(':',1)[0]
            return subset=='all' or k in sig or (unstranded in sig and len(by_unstranded[unstranded])==1)
          lengths=[v[0] for k,v in circ.items() if v[1]>=mincirc and selected(k)]
          exon_counts=[count_exons(k) for k,v in circ.items() if v[1]>=mincirc and selected(k)]
          row[subset+'_n_circs']=len(lengths);row[subset+'_median_span_bp']=statistics.median(lengths) if lengths else '';row[subset+'_median_exon_count']=statistics.median(exon_counts) if exon_counts else ''
      rows.append(row)
fields=['dataset','tissue','original_tissue','sample','junction_status']+[k+'_reads' for k in events]+['unc13a_junction_depth','unc13a_cryptic_fraction','stmn2_junction_depth','stmn2_cryptic_fraction','circ_status','all_n_circs','all_median_span_bp','all_median_exon_count','cqtl_n_circs','cqtl_median_span_bp','cqtl_median_exon_count','circ_eqtl_n_circs','circ_eqtl_median_span_bp','circ_eqtl_median_exon_count']
write('sample_circ_splicing.tsv',rows,fields)
for dataset,tissue in sorted({(r['dataset'],r['tissue']) for r in rows}):
  tissue_rows=[r for r in rows if r['dataset']==dataset and r['tissue']==tissue]
  # Three views: cryptic splicing, sample-median span, sample-median exon count.
  fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
  for ax,gene in zip(axes,['unc13a','stmn2']):
    fractions=[float(r[gene+'_cryptic_fraction']) for r in tissue_rows if r.get(gene+'_cryptic_fraction','')!='']
    if fractions:
      ax.hist([v*100 for v in fractions],bins=[0,1,5,10,25,50,75,100.00001],color='#426d89',edgecolor='white')
      ax.set_xlabel('Cryptic junction fraction (%)');ax.set_ylabel('Samples'); ax.set_title(gene.upper()+f' | n={len(fractions)}')
    else:
      ax.axis('off');ax.text(.5,.5,'No samples with sufficient junction coverage',ha='center',va='center',wrap=True,transform=ax.transAxes);ax.set_title(gene.upper()+' | n=0')
  fig.suptitle(dataset+' | '+tissue+' — cryptic splicing');fig.savefig(out/f'{dataset}_{tissue}_splicing_histogram.png',dpi=170);plt.close(fig)
  for metric,label,tag in [('median_span_bp','Median circRNA genomic span (bp)','span'),('median_exon_count','Median annotated exons per circRNA','exons')]:
    fig,axes=plt.subplots(2,3,figsize=(16,7),constrained_layout=True)
    for i,gene in enumerate(['unc13a','stmn2']):
      for j,subset in enumerate(['all','cqtl','circ_eqtl']):
        matched=[r for r in tissue_rows if r.get(subset+'_'+metric,'')!='' and r.get(gene+'_cryptic_fraction','')!='']
        pairs=[(float(r[subset+'_'+metric]),float(r[gene+'_cryptic_fraction'])) for r in matched]
        n=len(pairs);rho=pval=''
        if n>=3 and len({x for x,y in pairs})>1 and len({y for x,y in pairs})>1:
          res=spearmanr([x for x,y in pairs],[y for x,y in pairs]);rho=float(res.statistic);pval=float(res.pvalue)
        stats.append(dict(dataset=dataset,tissue=tissue,gene=gene,subset=subset,metric=metric,n=n,spearman_rho=rho,p_value=pval))
        ax=axes[i,j]; title=gene.upper()+' | '+{'all':'All detected','cqtl':'cQTL','circ_eqtl':'circ eQTL'}[subset]+f' | n={n}'
        ax.set_title(title)
        if not pairs:
          ax.axis('off');ax.text(.5,.5,'No matched samples with sufficient coverage',ha='center',va='center',wrap=True,transform=ax.transAxes);continue
        zero=[x for x,y in pairs if y==0];positive=[x for x,y in pairs if y>0]
        vals=[x for x,y in pairs]; lo=min(vals);hi=max(vals)
        if lo==hi: lo-=.5;hi+=.5
        import numpy as np
        bins=np.linspace(lo,hi,min(11,max(3,int(math.sqrt(n))+1)))
        ax.hist([zero,positive],bins=bins,stacked=True,color=['#a9b9c5','#d78758'],label=[f'No cryptic reads (n={len(zero)})',f'Cryptic reads detected (n={len(positive)})'],edgecolor='white')
        cryptic_counts=[(int(r['unc13a_128_reads'])+int(r['unc13a_178_reads']) if gene=='unc13a' else int(r['stmn2_cryptic_reads'])) for r in matched]
        detected_counts=[v for v in cryptic_counts if v>0]
        mean_all=statistics.mean(cryptic_counts)
        mean_detected=f'{statistics.mean(detected_counts):.2f}' if detected_counts else 'NA'
        ax.text(.98,.77,f'Mean cryptic reads: {mean_all:.2f} per sample\nAmong detected: {mean_detected}',transform=ax.transAxes,ha='right',va='top',fontsize=7,color='#444444',bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
        ax.set_xlabel(label);ax.set_ylabel('Samples');ax.legend(fontsize=7,frameon=False)
    fig.suptitle(dataset+' | '+tissue+' — sample distributions');fig.savefig(out/f'{dataset}_{tissue}_{tag}_histogram.png',dpi=170);plt.close(fig)
write('correlations.tsv',stats,['dataset','tissue','gene','subset','metric','n','spearman_rho','p_value'])
with (out/'issues.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['file','issue']);w.writerows(issues)
print('Output:',out,'samples:',len(rows),'issues:',len(issues),'minimum circ reads:',mincirc,'minimum competing junction reads:',minj)
PYTHON
