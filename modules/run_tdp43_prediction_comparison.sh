#!/bin/bash
#SBATCH --job-name=tdp43_prediction_comparison
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/tdp43_prediction.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/tdp43_prediction.err
set -euo pipefail
module purge
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1
module load deepTools
export INPUT_DIR="${INPUT_DIR:-$HOME/donglab/data/target_ALS/QTL/tdp43_cryptic_splicing_events}"
export OUTDIR="${OUTDIR:-$INPUT_DIR/prediction_comparison}"
export N_FOLDS="${N_FOLDS:-5}"
export N_BOOT="${N_BOOT:-2000}"
export SEED="${SEED:-1234}"
export MIN_DEPTH="${MIN_DEPTH:-10}"
mkdir -p "$OUTDIR"
python - <<'PY'
import os, csv, math, random, statistics
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import rankdata, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source=Path(os.environ['INPUT_DIR'])/'tdp43_pathology_sample_matches.tsv'
out=Path(os.environ['OUTDIR']); genes=('unc13a','stmn2','synj2','phf2')
cats=('Absent','Sparse','Moderate','Frequent')
fold_max=int(os.environ['N_FOLDS']); boot_n=int(os.environ['N_BOOT']); seed=int(os.environ['SEED']); min_depth=int(os.environ['MIN_DEPTH'])
if fold_max<2 or boot_n<1 or min_depth<1: raise ValueError('N_FOLDS >= 2, N_BOOT >= 1, MIN_DEPTH >= 1 required')
with source.open(newline='',encoding='utf-8-sig') as f:
 reader=csv.DictReader(f,delimiter='\t'); data=list(reader); columns=set(reader.fieldnames or [])
required={'dataset','tissue','original_tissue','sample','subject_id','pathology_match','junction_status','tdp43_neuronal','tdp43_glial'}
required |= {g+'_junction_depth' for g in genes} | {'unc13a_128_reads','unc13a_178_reads','unc13a_skip_reads'} | {g+'_'+suffix+'_reads' for g in genes[1:] for suffix in ('cryptic','canonical')}
if not required<=columns: raise ValueError('Input missing columns: '+', '.join(sorted(required-columns)))
def save(name,rows,fields):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def fraction(r,g):
 try:
  if r['junction_status']!='ok': return None
  depth=int(r[g+'_junction_depth'])
  if depth<min_depth: return None
  c=(int(r['unc13a_128_reads'])+int(r['unc13a_178_reads'])) if g=='unc13a' else int(r[g+'_cryptic_reads'])
  e=int(r['unc13a_skip_reads']) if g=='unc13a' else int(r[g+'_canonical_reads'])
  if c<0 or e<0 or c+e!=depth: raise ValueError('inconsistent junction counts')
  return c/depth
 except (ValueError,TypeError,KeyError) as e:
  raise ValueError('Invalid counts for %s %s %s: %s'%(r.get('sample'),r.get('tissue'),g,e))
def auc(y,p):
 y=np.asarray(y); p=np.asarray(p); a=y==1; n1=int(a.sum());n0=len(y)-n1
 if not n1 or not n0: return float('nan')
 return float((rankdata(p)[a].sum()-n1*(n1+1)/2)/(n1*n0))
def fit_logistic(X,y):
 # Fixed ridge strength, training-set scaling, no penalty on intercept.
 X=np.asarray(X,float); y=np.asarray(y,float); mu=X.mean(axis=0); sd=X.std(axis=0); sd[sd==0]=1
 X=(X-mu)/sd; D=np.column_stack((np.ones(len(y)),X)); lam=1.0
 def objective(b):
  eta=D@b
  return np.logaddexp(0,eta).sum()-y@eta+0.5*lam*np.dot(b[1:],b[1:])
 def gradient(b):
  grad=D.T@(expit(D@b)-y);grad[1:]+=lam*b[1:];return grad
 result=minimize(objective,np.zeros(D.shape[1]),jac=gradient,method='BFGS')
 if not result.success and np.linalg.norm(result.jac)>1e-4: raise RuntimeError('Logistic optimization failed: '+str(result.message))
 return mu,sd,result.x
def predict(model,X):
 mu,sd,b=model;return expit(b[0]+((np.asarray(X)-mu)/sd)@b[1:])
def folds(y,k,rng):
 assignments=np.empty(len(y),dtype=int)
 for label in (0,1):
  ids=np.flatnonzero(y==label).tolist();rng.shuffle(ids)
  for j,idx in enumerate(ids): assignments[idx]=j%k
 return assignments

audit=[]; feasibility=[]; metrics=[]; predictions=[]; bootstrap=[]; statuses=[]; severity=[]; ordinal=[]
tissues=sorted({r['tissue'] for r in data if r['dataset']=='target_ALS'})
for tissue in tissues:
 for cell in ('neuronal','glial'):
  records=[r for r in data if r['dataset']=='target_ALS' and r['tissue']==tissue and r['subject_id']]
  bysubject=defaultdict(list)
  for r in records: bysubject[r['subject_id']].append(r)
  chosen=[]
  for subject,rr in sorted(bysubject.items()):
   labels={r['tdp43_'+cell] for r in rr if r['tdp43_'+cell] in cats}
   if len(labels)!=1:
    audit.append(dict(tissue=tissue,cell_type=cell,subject_id=subject,status='missing_score' if not labels else 'conflicting_scores',candidate_samples=len(rr)));continue
   label=next(iter(labels)); candidates=[r for r in rr if r['tdp43_'+cell]==label]
   # Deterministic selection independent of gene-level junction coverage.
   # Prefer exact pathology specimen; otherwise select lexicographically.
   r=sorted(candidates,key=lambda a:(a['pathology_match']!='exact_sample',a['sample'],a['original_tissue']))[0]
   v=[fraction(r,g) for g in genes]
   chosen.append(dict(subject=subject,sample=r['sample'],original_tissue=r['original_tissue'],pathology_match=r['pathology_match'],category=label,y=int(label!='Absent'),values=v))
   audit.append(dict(tissue=tissue,cell_type=cell,subject_id=subject,sample=r['sample'],original_tissue=r['original_tissue'],category=label,pathology_match=r['pathology_match'],status='selected',candidate_samples=len(candidates),stmn2_available=int(v[1] is not None),four_available=int(all(x is not None for x in v))))
  for label in cats:
   items=[r for r in chosen if r['category']==label]
   feasibility.append(dict(tissue=tissue,cell_type=cell,category=label,selected_subjects=len(items),stmn2_available=sum(r['values'][1] is not None for r in items),all_four_available=sum(all(x is not None for x in r['values']) for r in items)))
  eligible=[r for r in chosen if all(x is not None for x in r['values'])]
  y=np.array([r['y'] for r in eligible],dtype=int)
  counts=Counter(y); k=min(fold_max,counts.get(0,0),counts.get(1,0)); key=dict(tissue=tissue,cell_type=cell,n_subjects=len(eligible),n_absent=counts.get(0,0),n_positive=counts.get(1,0))
  if k<2:
   statuses.append(dict(**key,status='insufficient_both_classes',detail='Need at least two subjects in each binary class with all four events measurable.'));continue
  X=np.asarray([r['values'] for r in eligible],float); assignment=folds(y,k,random.Random(seed)); pred=np.full((len(y),4),np.nan); fold_error=None; composite_error=[]
  for fold in range(k):
   train=assignment!=fold; test=~train
   # Composite parameters estimated only from absent training subjects, consistent
   # with the user's current Absent-reference score; all four events required.
   ref=X[train & (y==0)]; mu=ref.mean(axis=0); sd=ref.std(axis=0,ddof=1) if len(ref)>1 else np.full(4,np.nan)
   if np.any(~np.isfinite(sd)) or np.any(sd<=0):
    bad=','.join(genes[i] for i in range(4) if not np.isfinite(sd[i]) or sd[i]<=0)
    composite_error.append('Fold %d absent training reference has fewer than two values or zero SD: %s'%(fold+1,bad))
   try:
    pred[test,0]=predict(fit_logistic(X[train,1:2],y[train]),X[test,1:2])
    if not np.any(~np.isfinite(sd)) and np.all(sd>0): pred[test,1]=((X[test]-mu)/sd).mean(axis=1)
    pred[test,2]=predict(fit_logistic(X[train],y[train]),X[test])
    pred[test,3]=X[test].mean(axis=1)
   except Exception as e: fold_error='Fold %d: %s'%(fold+1,e);break
  if fold_error:
   statuses.append(dict(**key,status='invalid_training_fold',detail=fold_error));continue
  if composite_error: pred[:,1]=np.nan
  statuses.append(dict(**key,status='ok_composite_unavailable' if composite_error else 'ok',detail=('; '.join(composite_error) if composite_error else f'{k}-fold stratified subject-level cross-validation; complete cases')))
  names=('stmn2','equal_weight_composite','four_event_model','equal_weight_fraction')
  for i,r in enumerate(eligible):
   predictions.append(dict(tissue=tissue,cell_type=cell,subject_id=r['subject'],sample=r['sample'],original_tissue=r['original_tissue'],pathology_match=r['pathology_match'],category=r['category'],binary_outcome=r['y'],fold=assignment[i]+1,stmn2_fraction=X[i,1],stmn2_prediction=pred[i,0],composite_score=pred[i,1],four_event_prediction=pred[i,2],equal_weight_fraction=pred[i,3]))
  rng=np.random.default_rng(seed); distributions={n:[] for n in names}; diffs={'stmn2_minus_composite':[],'four_event_minus_stmn2':[],'stmn2_minus_equal_fraction':[]}
  # Stratified subject bootstrap resamples the already held-out predictions;
  # it quantifies uncertainty conditional on these fitted CV predictions.
  neg=np.flatnonzero(y==0);pos=np.flatnonzero(y==1)
  for b in range(boot_n):
   ix=np.concatenate((rng.choice(neg,len(neg),replace=True),rng.choice(pos,len(pos),replace=True)))
   vals=[auc(y[ix],pred[ix,j]) for j in range(4)]
   for name,v in zip(names,vals): distributions[name].append(v)
   diffs['stmn2_minus_composite'].append(vals[0]-vals[1]);diffs['four_event_minus_stmn2'].append(vals[2]-vals[0]);diffs['stmn2_minus_equal_fraction'].append(vals[0]-vals[3])
  for j,name in enumerate(names):
   if not np.all(np.isfinite(pred[:,j])): continue
   lo,hi=np.percentile(distributions[name],[2.5,97.5]);metrics.append(dict(**key,model=name,auc=auc(y,pred[:,j]),ci_low=lo,ci_high=hi,folds=k,bootstrap_replicates=boot_n))
   ordered=np.array([cats.index(r['category']) for r in eligible])
   rho=float(spearmanr(ordered,pred[:,j]).statistic) if len(set(ordered))>1 and len(set(pred[:,j]))>1 else float('nan')
   ordinal.append(dict(tissue=tissue,cell_type=cell,model=name,n=len(eligible),spearman_rho=rho,interpretation='descriptive association with ordered categories; not ordinal model validation'))
   for cat in cats:
    v=[pred[i,j] for i,r in enumerate(eligible) if r['category']==cat]
    severity.append(dict(tissue=tissue,cell_type=cell,category=cat,model=name,n=len(v),median_prediction=statistics.median(v) if v else ''))
  for name,v in diffs.items():
   if not np.all(np.isfinite(v)): continue
   estimate={'stmn2_minus_composite':lambda:auc(y,pred[:,0])-auc(y,pred[:,1]),'four_event_minus_stmn2':lambda:auc(y,pred[:,2])-auc(y,pred[:,0]),'stmn2_minus_equal_fraction':lambda:auc(y,pred[:,0])-auc(y,pred[:,3])}[name]();lo,hi=np.percentile(v,[2.5,97.5])
   bootstrap.append(dict(**key,comparison=name,auc_difference=estimate,ci_low=lo,ci_high=hi,bootstrap_replicates=boot_n))
  fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
  for j,name in enumerate(names):
   if not np.all(np.isfinite(pred[:,j])): continue
   ps=pred[:,j]; thresholds=np.r_[np.inf,np.sort(np.unique(ps))[::-1],-np.inf];fpr=[np.mean(ps[y==0]>=t) for t in thresholds];tpr=[np.mean(ps[y==1]>=t) for t in thresholds]
   axes[0].plot(fpr,tpr,color={'stmn2':'#0072B2','equal_weight_composite':'#D55E00','four_event_model':'#009E73','equal_weight_fraction':'#CC79A7'}[name],label=f"{dict(stmn2='STMN2 alone', equal_weight_composite='Equal-weight Z score', four_event_model='Fitted four-event model', equal_weight_fraction='Equal-weight raw fractions')[name]}: AUC {auc(y,ps):.2f}")
  axes[0].plot([0,1],[0,1],color='grey',ls='--');axes[0].set(xlabel='False positive rate',ylabel='True positive rate',title='Held-out predictions; same subjects');axes[0].legend(fontsize=8)
  for i,cat in enumerate(cats):
   v=[pred[j,0] for j,r in enumerate(eligible) if r['category']==cat]
   if v:
    axes[1].boxplot([v],positions=[i],widths=.5,showfliers=False,patch_artist=True,boxprops=dict(facecolor='#E6C68A',color='#504334'),medianprops=dict(color='#504334'),whiskerprops=dict(color='#504334'),capprops=dict(color='#504334'));axes[1].scatter([i]*len(v),v,s=12,alpha=.55,color='#765B3D')
  axes[1].set_xticks(range(4));axes[1].set_xticklabels([cat+'\n(n='+str(sum(r['category']==cat for r in eligible))+')' for cat in cats],rotation=0);axes[1].set(xlabel=cell.title()+' pTDP-43 category',ylabel='Held-out STMN2 predicted probability',title='Original ordered categories')
  fig.suptitle(f'{tissue} | {cell} | Absent vs any inclusions');fig.savefig(out/f'{tissue}_{cell}_prediction_comparison.png',dpi=160);plt.close(fig)

save('subject_selection_audit.tsv',audit,['tissue','cell_type','subject_id','sample','original_tissue','category','pathology_match','status','candidate_samples','stmn2_available','four_available'])
save('feasibility.tsv',feasibility,['tissue','cell_type','category','selected_subjects','stmn2_available','all_four_available'])
common=['tissue','cell_type','n_subjects','n_absent','n_positive']
save('analysis_status.tsv',statuses,common+['status','detail'])
save('auc_metrics.tsv',metrics,common+['model','auc','ci_low','ci_high','folds','bootstrap_replicates'])
save('auc_paired_differences.tsv',bootstrap,common+['comparison','auc_difference','ci_low','ci_high','bootstrap_replicates'])
save('held_out_predictions.tsv',predictions,['tissue','cell_type','subject_id','sample','original_tissue','pathology_match','category','binary_outcome','fold','stmn2_fraction','stmn2_prediction','composite_score','four_event_prediction','equal_weight_fraction'])
save('ordered_category_summary.tsv',severity,['tissue','cell_type','category','model','n','median_prediction'])
save('ordinal_association.tsv',ordinal,['tissue','cell_type','model','n','spearman_rho','interpretation'])
print('Results:',out,'| completed analyses:',sum(r['status'].startswith('ok') for r in statuses),'| attempted:',len(statuses))
PY
