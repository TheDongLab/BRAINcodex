#!/bin/bash
#SBATCH --job-name=tdp43_cryptic_splicing
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/tdp43_cryptic_splicing_events.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/tdp43_cryptic_splicing_events.err

set -euo pipefail
module purge

module load R

unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
module load deepTools
export TARGET_ROOT="${TARGET_ROOT:-$HOME/donglab/data/target_ALS}"
export ANSWER_ROOT="${ANSWER_ROOT:-$HOME/donglab/data/answer_ALS}"
export OUTDIR="${OUTDIR:-$TARGET_ROOT/QTL/tdp43_cryptic_splicing_events}"
export GTF="${GTF:-$HOME/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf}"
export TDP_PATHOLOGY="${TDP_PATHOLOGY:-$TARGET_ROOT/collections.postmortem_tissue_core.semiquantitative_tdp43_data.csv}"
export TARGET_METADATA="${TARGET_METADATA:-$TARGET_ROOT/targetALS_rnaseq_metadata.csv}"
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
events={'unc13a_128':('chr19',17642541,17642844,'-'),'unc13a_178':('chr19',17642591,17642844,'-'),'unc13a_skip':('chr19',17641556,17642844,'-'),'stmn2_cryptic':('chr8',79611214,79616821,'+'),'stmn2_canonical':('chr8',79611214,79636801,'+'),'synj2_cryptic':('chr6',158017290,158019983,'+'),'synj2_canonical':('chr6',158017290,158028755,'+'),'phf2_cryptic':('chr9',93660560,93661347,'+'),'phf2_canonical':('chr9',93660560,93662906,'+')}
def write(name,rows,fields):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
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
    for sd in sorted(x for x in processed.iterdir() if x.is_dir()):
      row={'dataset':dataset,'tissue':tissue,'original_tissue':original_tissue,'sample':sd.name}; junctions=list((sd/'leafcutter'/'juncs').glob('*.junc')); counts=Counter(); jstatus='ok'
      if len(junctions)!=1: jstatus='missing' if not junctions else 'multiple_files'
      else:
        try:
          with junctions[0].open() as f:
            for line in f:
              a=line.split()
              if len(a)<6: continue
              chrom='chr'+re.sub(r'^chr','',a[0],flags=re.I)
              if chrom not in ('chr6','chr8','chr9','chr19'): continue
              try: key=(chrom,int(a[1]),int(a[2]),a[5]); n=int(a[4])
              except ValueError: continue
              for name,coord in events.items():
                if key==coord: counts[name]+=n
        except OSError as e: jstatus='read_error';issues.append([str(junctions[0]),str(e)])
      row['junction_status']=jstatus
      for name in events: row[name+'_reads']=counts[name] if jstatus=='ok' else ''
      if jstatus=='ok':
        for gene,crypt,canonical in [('unc13a',counts['unc13a_128']+counts['unc13a_178'],counts['unc13a_skip']),('stmn2',counts['stmn2_cryptic'],counts['stmn2_canonical']),('synj2',counts['synj2_cryptic'],counts['synj2_canonical']),('phf2',counts['phf2_cryptic'],counts['phf2_canonical'])]:
          depth=crypt+canonical;row[gene+'_junction_depth']=depth;row[gene+'_cryptic_fraction']=crypt/depth if depth>=minj and depth>0 else ''
      rows.append(row)
fields=['dataset','tissue','original_tissue','sample','junction_status']+[k+'_reads' for k in events]+[item for gene in ['unc13a','stmn2','synj2','phf2'] for item in (gene+'_junction_depth',gene+'_cryptic_fraction')]
write('sample_cryptic_splicing.tsv',rows,fields)
for dataset,tissue in sorted({(r['dataset'],r['tissue']) for r in rows}):
  tissue_rows=[r for r in rows if r['dataset']==dataset and r['tissue']==tissue]
  # Three views: cryptic splicing, sample-median span, sample-median exon count.
  fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
  for ax,gene in zip(axes.flat,['unc13a','stmn2','synj2','phf2']):
    fractions=[float(r[gene+'_cryptic_fraction']) for r in tissue_rows if r.get(gene+'_cryptic_fraction','')!='']
    if fractions:
      # Zero is a discrete category; positive fractions use explicit unequal-width bins.
      bins=[0,0.1,0.5,1,5,10,25,50,75,100]
      labels=['0','>0–0.1','>0.1–0.5','>0.5–1','>1–5','>5–10','>10–25','>25–50','>50–75','>75–100']
      counts=[sum(v==0 for v in fractions)]
      percentages=[v*100 for v in fractions if v>0]
      counts.extend(sum(lo < v <= hi for v in percentages) for lo,hi in zip(bins[:-1],bins[1:]))
      ax.bar(range(len(labels)),counts,color=['#c44e52']+['#426d89']*(len(labels)-1),edgecolor='white')
      ax.set_xticks(range(len(labels)));ax.set_xticklabels(labels,rotation=45,ha='right',fontsize=8)
      ax.set_xlabel('Cryptic junction fraction (%)');ax.set_ylabel('Samples'); ax.set_title(gene.upper()+f' | n={len(fractions)}')
    else:
      ax.axis('off');ax.text(.5,.5,'No samples with sufficient junction coverage',ha='center',va='center',wrap=True,transform=ax.transAxes);ax.set_title(gene.upper()+' | n=0')
  fig.suptitle(dataset+' | '+tissue+' — cryptic splicing');fig.savefig(out/f'{dataset}_{tissue}_splicing_histogram.png',dpi=170);plt.close(fig)
# Link RNA samples to TDP-43 pathology using the RNA metadata subject ID and the same tissue region. An exact pathology sample ID is preferred. Subject + region matches may represent different specimens and are labelled accordingly.
def canonical_pathology_tissue(value):
  text=re.sub(r'[^a-z0-9]+',' ',str(value).lower()).strip()
  if 'cerebell' in text: return 'Cerebellum'
  if 'frontal' in text: return 'Frontal_Cortex'
  if 'motor' in text and ('cortex' in text or 'ba4' in text): return 'Motor_Cortex'
  if ('spinal' in text or 'cord' in text) and 'cervical' in text: return 'Cervical_Spinal_Cord'
  if ('spinal' in text or 'cord' in text) and ('lumbar' in text or 'lumbosacral' in text): return 'Lumbar_Spinal_Cord'
  return ''
def normid(value): return str(value or '').strip().upper()
def normsample(value): return normid(value).replace('_','-')
def read_csv_checked(path):
  with path.open(newline='',encoding='utf-8-sig') as f:
    reader=csv.DictReader(f)
    return list(reader),reader.fieldnames or []
pathology_path=Path(os.environ['TDP_PATHOLOGY']); metadata_path=Path(os.environ['TARGET_METADATA'])
path_index=defaultdict(list); meta_index=defaultdict(set)
if not pathology_path.is_file(): issues.append([str(pathology_path),'TDP-43 pathology CSV missing'])
else:
  pathology,columns=read_csv_checked(pathology_path)
  required={'Subject ID','Sample ID','Tissue Source','P Tdp 43 Inclusions Neuronal','P Tdp 43 Inclusions Glial'}
  if not required.issubset(columns): issues.append([str(pathology_path),'missing pathology columns: '+', '.join(sorted(required-set(columns)))])
  else:
    for record in pathology:
      subject=normid(record['Subject ID']); tissue=canonical_pathology_tissue(record['Tissue Source'])
      if subject and tissue: path_index[(subject,tissue)].append(record)
if not metadata_path.is_file(): issues.append([str(metadata_path),'Target ALS RNA metadata missing; pathology matching unavailable'])
else:
  meta,columns=read_csv_checked(metadata_path)
  lower={k.strip().lower():k for k in columns}
  samplecol=lower.get('externalsampleid'); subjectcol=lower.get('externalsubjectid')
  if not samplecol or not subjectcol: issues.append([str(metadata_path),'missing externalsampleid or externalsubjectid; available columns: '+', '.join(columns)])
  else:
    for record in meta:
      sample=normid(record.get(samplecol)); subject=normid(record.get(subjectcol))
      if sample and subject: meta_index[normsample(sample)].add((sample,subject))
missing_labels={'','N/D','NA','N/A','NOT AVAILABLE','UNKNOWN','NOT DONE'}
def score(value):
  value=str(value or '').strip()
  if value.upper() in missing_labels: return ''
  allowed={'absent':'Absent','sparse':'Sparse','moderate':'Moderate','frequent':'Frequent'}
  if value.lower() not in allowed:
    issues.append([str(pathology_path),'Unrecognized inclusion score: '+value]); return ''
  return allowed[value.lower()]
linked=[]
for row in rows:
  if row['dataset']!='target_ALS': continue
  sample=normid(row['sample']); entries=meta_index.get(normsample(sample),set()); status=''; record=None
  subjects={subject for metadata_sample,subject in entries}; metadata_ids={metadata_sample for metadata_sample,subject in entries}
  if not entries: status='rna_sample_not_in_metadata'
  elif len(subjects)>1: status='ambiguous_rna_subject'
  elif len(metadata_ids)>1: status='ambiguous_normalized_rna_sample'
  else:
    subject=next(iter(subjects)); candidates=path_index.get((subject,row['tissue']),[])
    exact=[r for r in candidates if normsample(r['Sample ID'])==normsample(sample)]
    if len(exact)==1: record=exact[0];status='exact_sample'
    elif len(exact)>1: status='ambiguous_exact_sample'
    elif len(candidates)==1: record=candidates[0];status='subject_tissue_different_specimen'
    elif len(candidates)>1: status='ambiguous_subject_tissue'
    else: status='no_pathology_for_subject_tissue'
  result=dict(row);result.update(metadata_sample_id=next(iter(metadata_ids)) if len(metadata_ids)==1 else '',subject_id=next(iter(subjects)) if len(subjects)==1 else '',pathology_match=status,pathology_sample_id='',pathology_tissue_source='',tdp43_neuronal='',tdp43_glial='')
  if record:
    result.update(pathology_sample_id=record['Sample ID'],pathology_tissue_source=record['Tissue Source'],tdp43_neuronal=score(record['P Tdp 43 Inclusions Neuronal']),tdp43_glial=score(record['P Tdp 43 Inclusions Glial']))
  linked.append(result)
linkfields=['dataset','tissue','original_tissue','sample','metadata_sample_id','subject_id','pathology_match','pathology_sample_id','pathology_tissue_source','tdp43_neuronal','tdp43_glial']+[x for x in fields if x not in ('dataset','tissue','original_tissue','sample')]
write('tdp43_pathology_sample_matches.tsv',linked,linkfields)
audit=Counter((r['tissue'],r['pathology_match']) for r in linked)
write('tdp43_match_summary.tsv',[dict(tissue=t,match_status=status,n=n) for (t,status),n in sorted(audit.items())],['tissue','match_status','n'])
# Select one sample per subject, tissue, gene and pathology-score type.
# Among eligible repeats, retain the sample with greatest competing-junction depth; ties use sample ID. This avoids counting a subject more than once per panel.
summary=[]; selected_samples=[]
markers=[('unc13a','UNC13A','#487895'),('stmn2','STMN2','#df8855'),('synj2','SYNJ2','#9d72af'),('phf2','PHF2','#c6a345')]
def event_reads(r,g):
 if g=='unc13a': return int(r['unc13a_128_reads'])+int(r['unc13a_178_reads'])
 return int(r[g+'_cryptic_reads'])
for tissue in sorted({r['tissue'] for r in linked}):
 for cell in ['neuronal','glial']:
  eligible=[r for r in linked if r['tissue']==tissue and r['tdp43_'+cell] and r['subject_id']]
  if not eligible: continue
  categories=['Absent','Sparse','Moderate','Frequent']; selected={}
  for gene,label,color in markers:
    candidates=defaultdict(list)
    for r in eligible:
      if r.get(gene+'_cryptic_fraction','')!='': candidates[r['subject_id']].append(r)
    selected[gene]=[]
    for subject,rr in sorted(candidates.items()):
      # Do not silently pick between conflicting pathology groups.
      if len({r['tdp43_'+cell] for r in rr})!=1:
        issues.append([subject,f'conflicting {cell} pathology scores for {tissue}; excluded for {gene}']);continue
      chosen=sorted(rr,key=lambda r:(-int(r[gene+'_junction_depth']),r['sample']))[0]
      selected[gene].append(chosen)
      selected_samples.append(dict(tissue=tissue,cell_type=cell,gene=gene,subject_id=subject,sample=chosen['sample'],category=chosen['tdp43_'+cell],junction_depth=chosen[gene+'_junction_depth'],pathology_match=chosen['pathology_match'],candidate_samples=len(rr)))
  junction_key=('Junction coordinates (.junc):\n'
    'UNC13A: cryptic chr19:17642541–17642844 + 17642591–17642844; canonical chr19:17641556–17642844\n'
    'STMN2: cryptic chr8:79611214–79616821; canonical chr8:79611214–79636801\n'
    'SYNJ2: cryptic chr6:158017290–158019983; canonical chr6:158017290–158028755\n'
    'PHF2: cryptic chr9:93660560–93661347; canonical chr9:93660560–93662906')
  def add_junction_key(fig): fig.text(.08,.015,junction_key,ha='left',va='bottom',fontsize=7,linespacing=1.3)
  fig,ax=plt.subplots(figsize=(12,8));fig.subplots_adjust(bottom=.43);width=.15
  for j,(gene,label,color) in enumerate(markers):
    for i,category in enumerate(categories):
      rr=[r for r in selected[gene] if r['tdp43_'+cell]==category];vals=[event_reads(r,gene) for r in rr];mean=statistics.mean(vals) if vals else None
      summary.append(dict(tissue=tissue,cell_type=cell,category=category,gene=gene,n=len(vals),mean_cryptic_reads=mean if mean is not None else '',sum_cryptic_reads=sum(vals) if vals else ''))
      if mean is not None:
        ax.bar(i+(j-1.5)*width,mean,width,color=color)
        ax.annotate(f'{mean:.2f}',(i+(j-1.5)*width,mean),xytext=(0,3),textcoords='offset points',ha='center',fontsize=7)
  from matplotlib.patches import Patch
  ax.set_xticks(range(4));ax.set_xticklabels(categories);ax.set_xlabel(cell.title()+' pTDP-43 inclusions');ax.set_ylabel('Mean event-supporting junction reads per subject');ax.set_title('Target ALS | '+tissue)
  handles=[Patch(color=color,label=label+'  n: '+', '.join(c+' '+str(sum(r['tdp43_'+cell]==c for r in selected[gene])) for c in categories)) for gene,label,color in markers]
  ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False,fontsize=8,ncol=1)
  add_junction_key(fig)
  fig.savefig(out/f'target_ALS_{tissue}_tdp43_{cell}_mean_cryptic_reads.png',dpi=170,bbox_inches='tight');plt.close(fig)
  fig,ax=plt.subplots(figsize=(12,8));fig.subplots_adjust(bottom=.43)
  for j,(gene,label,color) in enumerate(markers):
    for i,category in enumerate(categories):
      vals=[100*float(r[gene+'_cryptic_fraction']) for r in selected[gene] if r['tdp43_'+cell]==category]
      if vals: ax.boxplot([vals],positions=[i+(j-1.5)*width],widths=.12,patch_artist=True,boxprops=dict(facecolor=color,alpha=.65),medianprops=dict(color='black'),showfliers=True)
  ax.set_xticks(range(4));ax.set_xticklabels(categories);ax.set_xlabel(cell.title()+' pTDP-43 inclusions');ax.set_ylabel('Event junction fraction (%)');ax.set_title('Target ALS | '+tissue)
  ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False,fontsize=8,ncol=1)
  add_junction_key(fig)
  fig.savefig(out/f'target_ALS_{tissue}_tdp43_{cell}_combined_splicing.png',dpi=170,bbox_inches='tight');plt.close(fig)
  fig,ax=plt.subplots(figsize=(12,8));fig.subplots_adjust(bottom=.43)
  for j,(gene,label,color) in enumerate(markers):
    for i,category in enumerate(categories):
      rr=[r for r in selected[gene] if r['tdp43_'+cell]==category]; positive=sum(event_reads(r,gene)>0 for r in rr)
      ax.bar(i+(j-1.5)*width,len(rr)-positive,width,color=color,alpha=.25)
      ax.bar(i+(j-1.5)*width,positive,width,bottom=len(rr)-positive,color=color)
      if rr: ax.annotate(str(positive),(i+(j-1.5)*width,len(rr)),xytext=(0,3),textcoords='offset points',ha='center',fontsize=7)
  ax.set_xticks(range(4));ax.set_xticklabels(categories);ax.set_xlabel(cell.title()+' pTDP-43 inclusions');ax.set_ylabel('Subjects (solid = event detected)');ax.set_title('Target ALS | '+tissue)
  ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False,fontsize=8,ncol=1)
  add_junction_key(fig)
  fig.savefig(out/f'target_ALS_{tissue}_tdp43_{cell}_combined_detection.png',dpi=170,bbox_inches='tight');plt.close(fig)
write('tdp43_selected_subject_samples.tsv',selected_samples,['tissue','cell_type','gene','subject_id','sample','category','junction_depth','pathology_match','candidate_samples'])
write('tdp43_mean_cryptic_reads.tsv',summary,['tissue','cell_type','category','gene','n','mean_cryptic_reads','sum_cryptic_reads'])
with (out/'issues.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['file','issue']);w.writerows(issues)
print('Output:',out,'samples:',len(rows),'issues:',len(issues),'minimum circ reads:',mincirc,'minimum competing junction reads:',minj)
PYTHON

# Four neuronal pathology groups using the unmodified LeafViz plotting script.
PLOT_SCRIPT="${PLOT_SCRIPT:-$HOME/donglab/pipelines/scripts/QTL/plot_sQTL_leafviz_PSI.R}"
TDP_ALIGNED="${TDP_ALIGNED:-$TARGET_ROOT/aligned_rnaseq_tdp43_by_subject.csv}"
RSCRIPT="${RSCRIPT:-Rscript}"
PLOT_OUTDIR="${PLOT_OUTDIR:-$OUTDIR/leafviz_tdp43}"
mkdir -p "$PLOT_OUTDIR"
[[ -r "$PLOT_SCRIPT" ]] || { echo "Missing plotting script: $PLOT_SCRIPT" >&2; exit 1; }
[[ -r "$TDP_ALIGNED" ]] || { echo "Missing aligned pathology table: $TDP_ALIGNED" >&2; exit 1; }
command -v "$RSCRIPT" >/dev/null || { echo 'Rscript unavailable; load the appropriate R module' >&2; exit 1; }
# PSI-file anchor endpoints use .junc end + 1, as documented for UNC13A
# in the plotting script. Override these if your PSI files use another convention.
UNC13A_ANCHOR="${UNC13A_ANCHOR:-chr19:-:17641556-17642845}"
STMN2_ANCHOR="${STMN2_ANCHOR:-chr8:+:79611214-79636802}"
SYNJ2_ANCHOR="${SYNJ2_ANCHOR:-chr6:+:158017290-158028756}"
PHF2_ANCHOR="${PHF2_ANCHOR:-chr9:+:93660560-93662907}"
plot_one() {
  local tissue="$1" tissue_regex="$2" gene="$3" anchor="$4"
  local prefix="target_ALS_${tissue}_${gene}_tdp43_neuronal"
  echo "Plotting $gene in $tissue"
  "$RSCRIPT" "$PLOT_SCRIPT" --group-by tdp43 \
    --tissue "$tissue" --tissue-regex "$tissue_regex" \
    --anchor "$anchor" --gene "$gene" --metadata "$TARGET_METADATA" \
    --data-root "$TARGET_ROOT" --gtf "$GTF" --tdp43 "$TDP_ALIGNED" \
    --pathology "$TDP_PATHOLOGY" --outdir "$PLOT_OUTDIR" --prefix "$prefix"
}
# Default to the illustrated tissue; override PLOT_TISSUES with a space-separated list.
for tissue in ${PLOT_TISSUES:-Lumbar_Spinal_Cord}; do
  case "$tissue" in
    Lumbar_Spinal_Cord) tissue_regex='lumbar|lumbosacral' ;;
    Cervical_Spinal_Cord) tissue_regex='cervical' ;;
    Motor_Cortex) tissue_regex='motor|BA4' ;;
    Frontal_Cortex) tissue_regex='frontal' ;;
    Cerebellum) tissue_regex='cerebell' ;;
    *) echo "Unsupported plot tissue: $tissue" >&2; exit 1 ;;
  esac
  plot_one "$tissue" "$tissue_regex" UNC13A "$UNC13A_ANCHOR"
  plot_one "$tissue" "$tissue_regex" STMN2 "$STMN2_ANCHOR"
  plot_one "$tissue" "$tissue_regex" SYNJ2 "$SYNJ2_ANCHOR"
  plot_one "$tissue" "$tissue_regex" PHF2 "$PHF2_ANCHOR"
done
