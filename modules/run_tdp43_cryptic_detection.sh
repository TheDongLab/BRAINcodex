#!/bin/bash
#SBATCH --job-name=tdp43_ce_detection
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/tdp43_ce_detection.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/tdp43_ce_detection.err

set -euo pipefail
module purge
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
module load deepTools
export INPUT_DIR="${INPUT_DIR:-$HOME/donglab/data/target_ALS/QTL/tdp43_cryptic_splicing_events}"
export OUTDIR="${OUTDIR:-$INPUT_DIR/ce_detection}"
export MIN_JUNCTION_READS="${MIN_JUNCTION_READS:-10}"
mkdir -p "$OUTDIR"
python - <<'PYTHON'
import csv,os,math
from pathlib import Path
from collections import defaultdict,Counter
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source=Path(os.environ['INPUT_DIR'])/'tdp43_pathology_sample_matches.tsv'
out=Path(os.environ['OUTDIR']); minimum=int(os.environ['MIN_JUNCTION_READS'])
if minimum<1: raise ValueError('MIN_JUNCTION_READS must be positive')
genes=('unc13a','stmn2','synj2','phf2'); cats=('Absent','Sparse','Moderate','Frequent')
colors={'unc13a':'#487895','stmn2':'#df8855','synj2':'#9d72af','phf2':'#c6a345'}
with source.open(newline='',encoding='utf-8-sig') as f:
 reader=csv.DictReader(f,delimiter='\t');columns=set(reader.fieldnames or []);data=list(reader)
need={'dataset','tissue','original_tissue','sample','subject_id','junction_status','pathology_match','pathology_sample_id','tdp43_neuronal','tdp43_glial','unc13a_128_reads','unc13a_178_reads','unc13a_skip_reads'}|{g+'_junction_depth' for g in genes}|{g+'_'+x+'_reads' for g in genes[1:] for x in ('cryptic','canonical')}
if need-columns: raise ValueError('Missing input columns: '+', '.join(sorted(need-columns)))
def write(name,records,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(records)
def binom(k,n):
 if not n:return ('','','')
 return (k/n,0 if k==0 else beta.ppf(.025,k,n-k+1),1 if k==n else beta.ppf(.975,k+1,n-k))
def counts(r,g):
 if r['junction_status']!='ok':return None
 try:
  c=int(r['unc13a_128_reads'])+int(r['unc13a_178_reads']) if g=='unc13a' else int(r[g+'_cryptic_reads'])
  e=int(r['unc13a_skip_reads']) if g=='unc13a' else int(r[g+'_canonical_reads'])
  depth=int(r[g+'_junction_depth'])
  if c<0 or e<0 or c+e!=depth: raise ValueError('inconsistent counts')
  return c,e,depth
 except (KeyError,ValueError) as e: raise ValueError('Invalid junction counts for '+r['sample']+' '+g+': '+str(e))
audit=[];individual=[];category=[];binary=[]
for tissue in sorted({r['tissue'] for r in data if r['dataset']=='target_ALS'}):
 for cell in ('neuronal','glial'):
  bysubject=defaultdict(list)
  for r in data:
   if r['dataset']=='target_ALS' and r['tissue']==tissue and r['subject_id']: bysubject[r['subject_id']].append(r)
  chosen=[]
  for subject,rr in sorted(bysubject.items()):
   labels={r['tdp43_'+cell] for r in rr if r['tdp43_'+cell] in cats}
   if len(labels)!=1:
    audit.append(dict(tissue=tissue,cell_type=cell,subject_id=subject,status='missing_pathology' if not labels else 'conflicting_pathology',candidate_samples=len(rr)));continue
   cat=next(iter(labels)); candidates=[r for r in rr if r['tdp43_'+cell]==cat]
   # One sample chosen before examining any event counts.
   r=sorted(candidates,key=lambda x:(x['pathology_match']!='exact_sample',x['sample'],x['original_tissue']))[0]
   chosen.append(r);audit.append(dict(tissue=tissue,cell_type=cell,subject_id=subject,sample=r['sample'],original_tissue=r['original_tissue'],pathology_sample_id=r['pathology_sample_id'],pathology_match=r['pathology_match'],category=cat,status='selected',candidate_samples=len(candidates)))
  for gene in genes:
   obs=[]
   for r in chosen:
    count=counts(r,gene); c,e,depth=count if count else ('','','')
    evaluable=count is not None and depth>=minimum
    obs.append((r,c,e,depth,evaluable))
    individual.append(dict(tissue=tissue,cell_type=cell,gene=gene.upper(),subject_id=r['subject_id'],sample=r['sample'],original_tissue=r['original_tissue'],pathology_sample_id=r['pathology_sample_id'],pathology_match=r['pathology_match'],category=r['tdp43_'+cell],cryptic_reads=c,canonical_reads=e,depth=depth,evaluable=int(evaluable),detected_2_reads=int(c>=2) if evaluable else '',detected_1_read=int(c>=1) if evaluable else ''))
   for threshold in (2,1):
    for cat in cats:
     group=[x for x in obs if x[0]['tdp43_'+cell]==cat]; ok=[x for x in group if x[4]]; n=len(ok);k=sum(x[1]>=threshold for x in ok);rate,lo,hi=binom(k,n)
     category.append(dict(tissue=tissue,cell_type=cell,gene=gene.upper(),threshold_reads=threshold,category=cat,total_selected=len(group),evaluable=n,not_evaluable=len(group)-n,detected=k,not_detected=n-k,detection_rate=rate,ci_low=lo,ci_high=hi))
    neg=[x for x in obs if x[4] and x[0]['tdp43_'+cell]=='Absent'];pos=[x for x in obs if x[4] and x[0]['tdp43_'+cell] in cats[1:]]
    tp=sum(x[1]>=threshold for x in pos); fn=len(pos)-tp;fp=sum(x[1]>=threshold for x in neg);tn=len(neg)-fp
    sens,sl,sh=binom(tp,tp+fn);spec,cl,ch=binom(tn,tn+fp)
    binary.append(dict(tissue=tissue,cell_type=cell,gene=gene.upper(),threshold_reads=threshold,tp=tp,fn=fn,fp=fp,tn=tn,sensitivity=sens,sensitivity_ci_low=sl,sensitivity_ci_high=sh,specificity=spec,specificity_ci_low=cl,specificity_ci_high=ch,selected_subjects=len(chosen),not_evaluable=sum(not x[4] for x in obs)))
  # Plot prespecified two-read detection, not sensitivity-analysis threshold.
  fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True);width=.18
  for j,gene in enumerate(genes):
   for i,cat in enumerate(cats):
    row=next(r for r in category if r['tissue']==tissue and r['cell_type']==cell and r['gene']==gene.upper() and r['threshold_reads']==2 and r['category']==cat)
    if not row['evaluable']:continue
    rate=row['detection_rate'];x=i+(j-1.5)*width
    ax.bar(x,rate*100,width,color=colors[gene],label=gene.upper() if i==0 else None)
    ax.errorbar([x],[rate*100],yerr=[[100*(rate-row['ci_low'])],[100*(row['ci_high']-rate)]],fmt='none',ecolor='black',capsize=2,lw=.8)
    ax.annotate(str(row['detected'])+'/'+str(row['evaluable']),(x,rate*100),xytext=(0,3),textcoords='offset points',ha='center',fontsize=7)
  from matplotlib.patches import Patch
  ax.legend(handles=[Patch(color=colors[g],label=g.upper()) for g in genes],frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(.5,-.17))
  ax.set_xticks(range(4));ax.set_xticklabels(cats);ax.set_ylim(0,115);ax.set_ylabel('Subjects with detected cryptic junction (%)');ax.set_xlabel(cell.title()+' pTDP-43 inclusions');ax.set_title(tissue+' | '+cell+' | ≥2 cryptic reads; '+str(minimum)+' total junction reads')
  fig.savefig(out/f'{tissue}_{cell}_ce_detection.png',dpi=170,bbox_inches='tight');plt.close(fig)
# Coverage figures and matched-subject detection comparison.
# Event depth means cryptic + selected canonical junction reads, not library size.
paired=[]
for tissue,cell in sorted({(r['tissue'],r['cell_type']) for r in individual}):
 rows=[r for r in individual if r['tissue']==tissue and r['cell_type']==cell]
 by_subject=defaultdict(dict)
 for r in rows: by_subject[r['subject_id']][r['gene']]=r
 common=[v for v in by_subject.values() if all(g.upper() in v and v[g.upper()]['evaluable']==1 for g in genes)]
 fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
 for j,g in enumerate(genes):
  vals=[int(r['depth']) for r in rows if r['gene']==g.upper() and r['depth']!='']
  if vals:
   axes[0].boxplot([vals],positions=[j],widths=.5,showfliers=False,patch_artist=True,boxprops=dict(facecolor=colors[g],alpha=.6),medianprops=dict(color='black'))
   axes[0].scatter([j]*len(vals),vals,s=9,alpha=.3,color=colors[g])
  total=sum(r['gene']==g.upper() for r in rows); available=sum(r['gene']==g.upper() and r['evaluable']==1 for r in rows)
  axes[1].bar(j,100*available/total if total else 0,color=colors[g])
  axes[1].text(j,100*available/total+2 if total else 2,f'{available}/{total}',ha='center',fontsize=9)
 axes[0].set_yscale('symlog',linthresh=1);axes[0].axhline(minimum,color='black',ls='--',lw=1,label=f'Coverage threshold = {minimum}');axes[0].legend(fontsize=8)
 axes[0].set_xticks(range(4));axes[0].set_xticklabels([g.upper()+'\n(n='+str(sum(r['gene']==g.upper() and r['depth']!='' for r in rows))+')' for g in genes]);axes[0].set_ylabel('Selected junction reads (cryptic + canonical)');axes[0].set_title('Event coverage; all pathology categories')
 axes[1].set_xticks(range(4));axes[1].set_xticklabels([g.upper() for g in genes]);axes[1].set_ylim(0,115);axes[1].set_ylabel('Selected subjects evaluable (%)');axes[1].set_title('At least '+str(minimum)+' junction reads')
 fig.suptitle(tissue+' | '+cell+' — event coverage');fig.savefig(out/f'{tissue}_{cell}_event_coverage.png',dpi=170);plt.close(fig)
 fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True);width=.18
 for j,g in enumerate(genes):
  for i,cat in enumerate(cats):
   group=[v for v in common if v[g.upper()]['category']==cat]; n=len(group); k=sum(v[g.upper()]['detected_2_reads']==1 for v in group)
   rate,lo,hi=binom(k,n)
   paired.append(dict(tissue=tissue,cell_type=cell,category=cat,gene=g.upper(),n_common=n,detected=k,rate=rate,ci_low=lo,ci_high=hi))
   if n:
    x=i+(j-1.5)*width;ax.bar(x,rate*100,width,color=colors[g]);ax.errorbar([x],[rate*100],yerr=[[100*(rate-lo)],[100*(hi-rate)]],fmt='none',ecolor='black',capsize=2,lw=.8);ax.annotate(f'{k}/{n}',(x,rate*100),xytext=(0,3),textcoords='offset points',ha='center',fontsize=7)
 from matplotlib.patches import Patch
 ax.legend(handles=[Patch(color=colors[g],label=g.upper()) for g in genes],frameon=False,ncol=4,loc='upper center',bbox_to_anchor=(.5,-.17))
 ax.set_xticks(range(4));ax.set_xticklabels(cats);ax.set_ylim(0,115);ax.set_ylabel('Subjects with detected cryptic junction (%)');ax.set_xlabel(cell.title()+' pTDP-43 inclusions');ax.set_title(tissue+' | '+cell+' | same subjects evaluable for all four events')
 fig.savefig(out/f'{tissue}_{cell}_ce_detection_common_subjects.png',dpi=170,bbox_inches='tight');plt.close(fig)
write('common_subject_detection.tsv',paired,['tissue','cell_type','category','gene','n_common','detected','rate','ci_low','ci_high'])
write('subject_selection_audit.tsv',audit,['tissue','cell_type','subject_id','sample','original_tissue','pathology_sample_id','pathology_match','category','status','candidate_samples'])
write('subject_event_detection.tsv',individual,['tissue','cell_type','gene','subject_id','sample','original_tissue','pathology_sample_id','pathology_match','category','cryptic_reads','canonical_reads','depth','evaluable','detected_2_reads','detected_1_read'])
write('detection_by_category.tsv',category,['tissue','cell_type','gene','threshold_reads','category','total_selected','evaluable','not_evaluable','detected','not_detected','detection_rate','ci_low','ci_high'])
write('binary_performance.tsv',binary,['tissue','cell_type','gene','threshold_reads','tp','fn','fp','tn','sensitivity','sensitivity_ci_low','sensitivity_ci_high','specificity','specificity_ci_low','specificity_ci_high','selected_subjects','not_evaluable'])
print('Detection results:',out)
PYTHON
