#!/bin/bash
#SBATCH --job-name=rbp_junction_map
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/rbp_junction_map_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/rbp_junction_map_%j.err
set -euo pipefail
# Usage: sbatch run_rbp_junction_map.sh [--cohort both|target|answer] [--tissues TISSUE ...]
# Override Python with TASK_PYTHON; run --help for all options.
TASK_PYTHON=${TASK_PYTHON:-/home/zw529/donglab/pipelines/modules/miniconda3/bin/python}
if [[ ! -x "$TASK_PYTHON" ]]; then
    echo "ERROR: Python executable not found: $TASK_PYTHON" >&2
    exit 1
fi
"$TASK_PYTHON" -c 'import sys, numpy, pandas, scipy; print("Python:", sys.executable); print("numpy:", numpy.__version__, "pandas:", pandas.__version__, "scipy:", scipy.__version__)' || {
    echo "ERROR: Install numpy, pandas and scipy in this Python environment, or set TASK_PYTHON to a suitable Python executable." >&2
    exit 1
}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$TASK_PYTHON" -u - "$@" <<'PY'
import argparse, csv, gzip, json, math, re, sys, sqlite3, hashlib, os
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
try:
 from scipy.stats import t as tdist
 from scipy.linalg import qr
except ImportError:
 sys.exit('Missing scipy; set TASK_PYTHON to an environment containing numpy, pandas and scipy.')

HOME=Path.home(); BASE=HOME/'donglab/data'; REF=HOME/'donglab/references/genome/Homo_sapiens'
p=argparse.ArgumentParser(description='Cohort- and tissue-specific RBP expression / unannotated junction associations.')
p.add_argument('--force',action='store_true',help='Recompute completed tissues')
p.add_argument('--batch-size',type=int,default=512,help='RBP batch size')
p.add_argument('--cohort',choices=['both','target','answer'],default='both')
p.add_argument('--tissues',nargs='+',help='Exact tissue directory names; default: discover directories with RNAseq/Processed')
p.add_argument('--outdir',type=Path,default=BASE/'target_ALS/QTL/rbp_unannotated_junction_associations')
p.add_argument('--gtf',type=Path,default=REF/'UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf')
p.add_argument('--rbp-list',type=Path,default=REF/'RBP2GO/RNA_binding_proteins/RNABP_list.csv')
p.add_argument('--min-samples',type=int,default=20)
p.add_argument('--min-reads',type=int,default=5,help='Minimum cluster reads per sample for a valid usage measurement')
p.add_argument('--min-carriers',type=int,default=5,help='Minimum samples with at least one junction read')
p.add_argument('--min-tpm',type=float,default=1.0,help='Minimum median RBP TPM')
p.add_argument('--fdr',type=float,default=0.05)
p.add_argument('--max-map-edges',type=int,default=500,help='Maximum edges displayed in HTML; all significant results remain in tables')
a=p.parse_args()
if a.batch_size<1 or a.min_samples<5 or a.min_reads<1 or a.min_carriers<2 or not 0<a.fdr<=1: p.error('Invalid threshold')
for f in [a.gtf,a.rbp_list]:
 if not f.is_file(): sys.exit('Missing input: '+str(f))
def chrom(s): return re.sub(r'^chr','',str(s),flags=re.I)
def normid(s): return re.sub(r'[-_]','',str(s)).upper()
def op(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)
def bh(x):
 x=np.asarray(x); order=np.argsort(x); vals=x[order]*len(x)/np.arange(1,len(x)+1); vals=np.minimum.accumulate(vals[::-1])[::-1]; out=np.empty(len(x)); out[order]=np.minimum(vals,1); return out
# GTF exon pairs: PSI junction convention is upstream exon end and downstream exon start.
rbps=pd.read_csv(a.rbp_list,comment='#',usecols=['Gene_Name','Is_RBP'],dtype=str).fillna('')
rbps=rbps[rbps.Is_RBP.str.upper().str.strip()=='TRUE']
if 'Gene_Name' not in rbps: sys.exit('RBP list lacks Gene_Name')
symbols=set(rbps.Gene_Name.str.strip())-{''}
transcripts=defaultdict(list); gene_symbols=defaultdict(set)
with op(a.gtf) as f:
 for line in f:
  if line.startswith('#'): continue
  x=line.rstrip('\n').split('\t')
  if len(x)<9 or x[2]!='exon': continue
  attrs=dict(re.findall(r'(\S+)\s+"([^"]+)"',x[8])); gid=attrs.get('gene_id','').split('.')[0]; name=attrs.get('gene_name',''); tx=attrs.get('transcript_id','')
  if name in symbols and gid: gene_symbols[gid].add(name)
  if tx: transcripts[(chrom(x[0]),x[6],tx)].append((int(x[3]),int(x[4])))
annotated=set()
for (c,s,_),exons in transcripts.items():
 exons=sorted(set(exons))
 for left,right in zip(exons,exons[1:]):
  if left[1]<right[0]: annotated.add((c,s,left[1],right[0]))
del transcripts
print('RBP symbols:',len(symbols),'mapped Ensembl IDs:',len(gene_symbols),'annotated junctions:',len(annotated),flush=True)

def read_junc(path):
 d=Counter()
 with op(path) as f:
  for line in f:
   v=line.split()
   if len(v)<6 or v[5] not in ('+','-'): continue
   try: start,end,count=int(v[1]),int(v[2])+1,int(v[4])
   except ValueError: continue
   if count<0 or end<=start or end-start>100000: continue
   d[(chrom(v[0]),v[5],start,end)]+=count
 return d

def expression(path):
 df=pd.read_csv(path,sep='\t',usecols=['gene_id','TPM'],dtype={'gene_id':str})
 df['gene_id']=df.gene_id.str.split('.').str[0]; df['TPM']=pd.to_numeric(df.TPM,errors='coerce')
 df=df[df.gene_id.isin(gene_symbols)].groupby('gene_id',as_index=True).TPM.sum(min_count=1)
 return df.to_dict()

def loadmeta(cohort,root):
 if cohort=='target':
  m=pd.read_csv(root/'targetALS_rnaseq_metadata.csv',dtype=str,low_memory=False).fillna(''); key='externalsampleid'; subject='externalsubjectid'
  cov=['sex','age_at_death','rin','post_mortem_interval_in_hours','platform','subject_group']
 else:
  m=pd.read_csv(root/'aals_dataportal_datatable.csv',dtype=str,low_memory=False).fillna(''); key='Participant_ID'; subject='Participant_ID'
  cov=['SEX','Age_at_First_PBMC_Collection','Batch','SUBJECT_GROUP']
 if key not in m: raise ValueError('Metadata missing '+key)
 m['_key']=m[key].map(normid); m=m[m._key!=''].copy(); cols=[subject]+[c for c in cov if c in m]
 # Collapse only consistent fields; conflicting values remain missing and are audited.
 rows=[]
 for k,g in m.groupby('_key',sort=False):
  row={'_key':k}
  for c in cols:
   vals=sorted(set(v.strip() for v in g[c] if str(v).strip() and str(v).lower() not in ('nan','unknown','na','n/a')))
   row[c]=vals[0] if len(vals)==1 else ''
  rows.append(row)
 return pd.DataFrame(rows).set_index('_key'),subject,[c for c in cov if c in m]

def design(meta,covs):
 # Complete-case design; only covariates with sufficient nonmissing data are selected.
 keep=[]; encoded=[]
 for c in covs:
  z=meta[c].replace('',np.nan); z=z.str.strip().str.lower() if c.lower() in ('sex','subject_group') else z; numeric=pd.to_numeric(z,errors='coerce')
  numeric_col=c.lower() in ('age_at_death','rin','post_mortem_interval_in_hours','age_at_first_pbmc_collection')
  if numeric_col:
   if numeric.notna().sum()>=max(a.min_samples,int(.8*len(meta))) and numeric.nunique()>1: keep.append(c); encoded.append(numeric.rename(c))
  else:
   if z.notna().sum()>=max(a.min_samples,int(.8*len(meta))) and 1<z.nunique()<=10:
    keep.append(c); encoded.append(pd.get_dummies(z,prefix=c,drop_first=True,dtype=float).where(z.notna(),np.nan))
 if not encoded: return np.ones((len(meta),1)),np.ones(len(meta),dtype=bool),keep
 mat=pd.concat(encoded,axis=1); valid=mat.notna().all(axis=1).to_numpy(); X=np.column_stack([np.ones(len(meta)),mat.fillna(0).to_numpy(dtype=float)])
 return X,valid,keep

def run(cohort):
 root=BASE/('target_ALS' if cohort=='target' else 'answer_ALS'); metadata,subject,covs=loadmeta(cohort,root)
 dirs=sorted([x for x in root.iterdir() if x.is_dir() and (x/'RNAseq/Processed').is_dir()]); dirs=[x for x in dirs if not a.tissues or x.name in a.tissues]
 print(cohort,'tissues:',[x.name for x in dirs],flush=True)
 for tissue in dirs:
  out=a.outdir/cohort/tissue.name; out.mkdir(parents=True,exist_ok=True)
  marker=out/'COMPLETE.json'
  options={k:str(v) for k,v in vars(a).items() if k not in ('force','tissues','cohort')}
  fingerprint=hashlib.sha256(json.dumps(options,sort_keys=True).encode()).hexdigest()
  if marker.exists() and not a.force:
   old=json.loads(marker.read_text())
   if old.get('settings')==fingerprint:
    print('Skipping completed tissue:',cohort,tissue.name,flush=True); continue
   raise RuntimeError('Completed tissue has different settings; use --force: '+str(out))
  if marker.exists(): marker.unlink()
  audit=[]; samples=[]; seen=set(); processed=tissue/'RNAseq/Processed'
  for d in sorted(processed.iterdir()):
   if not d.is_dir(): continue
   key=normid(d.name); jf=d/'leafcutter/juncs'/f'{d.name}.junc'; ef=d/'normalization.tab'; reason=''
   if not jf.is_file() or not ef.is_file(): reason='missing_junction_or_expression'
   elif ef.stat().st_size==0 or sum(1 for _ in ef.open())<=1: reason='empty_expression_file'
   elif key not in metadata.index: reason='metadata_not_matched'
   elif not str(metadata.loc[key,subject]).strip(): reason='missing_or_conflicting_subject'
   else:
    sid=normid(metadata.loc[key,subject])
    if sid in seen: reason='duplicate_subject_excluded'
    else: seen.add(sid)
   audit.append({'sample':d.name,'metadata_key':key,'status':reason or 'included'})
   if not reason: samples.append((d.name,key,jf,ef))
  pd.DataFrame(audit,columns=['sample','metadata_key','status']).to_csv(out/'sample_audit.tsv',sep='\t',index=False)
  if len(samples)<a.min_samples: print(tissue.name,'too few matched samples',len(samples),flush=True); continue
  meta=metadata.loc[[s[1] for s in samples]].copy(); X,valid,covused=design(meta,covs)
  # Check rank after complete-case selection; avoid a false regression if design is singular.
  if valid.sum()<a.min_samples:
   print(tissue.name,'insufficient complete cases:',int(valid.sum()),flush=True); continue
  _,_,pivot=qr(X[valid],mode='economic',pivoting=True)
  rank=np.linalg.matrix_rank(X[valid]); selected=sorted(pivot[:rank]); dropped=X.shape[1]-rank
  if dropped: print('Dropping',dropped,'redundant model columns',flush=True)
  X=X[:,selected]
  print(cohort,tissue.name,'samples',len(samples),'complete',valid.sum(),'covariates',covused,flush=True)
  js=[]; ex=[]; pooled=Counter()
  for name,key,jf,ef in samples:
   j=read_junc(jf); js.append(j); pooled.update(j); ex.append(expression(ef))
  # Shared clusters are connected by shared splice sites, after pooled support >=5.
  groups=defaultdict(list)
  for (c,s,start,end),n in pooled.items():
   if n>=5: groups[(c,s)].append((start,end))
  clusters=[]
  for (c,s),intervals in groups.items():
   parent={}
   def find(x):
    parent.setdefault(x,x)
    if parent[x]!=x: parent[x]=find(parent[x])
    return parent[x]
   def union(x,y):
    x=find(x); y=find(y)
    if x!=y: parent[y]=x
   for start,stop in intervals: union(('L',start),('R',stop))
   sets=defaultdict(list)
   for start,stop in intervals: sets[find(('L',start))].append((start,stop))
   for ivs in sets.values(): clusters.append((c,s,ivs))
  print('Shared clusters:',len(clusters),'largest:',max((len(x[2]) for x in clusters),default=0),flush=True)
  geneids=sorted(gene_symbols); E=np.array([[row.get(g,np.nan) for g in geneids] for row in ex],dtype=float)
  genes=[i for i in range(len(geneids)) if np.isfinite(E[valid,i]).sum()>=a.min_samples and np.nanmedian(E[valid,i])>=a.min_tpm and np.nanstd(E[valid,i])>0]
  print('Eligible RBP genes:',len(genes),flush=True)
  dbpath=out/'working_associations.sqlite'
  if dbpath.exists(): dbpath.unlink()
  con=sqlite3.connect(str(dbpath)); con.execute('PRAGMA journal_mode=DELETE'); con.execute('PRAGMA synchronous=NORMAL'); con.execute('PRAGMA temp_store=FILE')
  cols=['rbp_ensembl','rbp_gene','junction_id','chrom','strand','start','end','n_samples','n_with_reads','partial_r','p_value','cluster_junctions']
  con.execute('CREATE TABLE pairs (rbp_ensembl TEXT, rbp_gene TEXT, junction_id TEXT, chrom TEXT, strand TEXT, start INTEGER, end INTEGER, n_samples INTEGER, n_with_reads INTEGER, partial_r REAL, p_value REAL, cluster_junctions INTEGER, fdr REAL)')
  tested_events=0; total_pairs=0
  logged=np.log2(E+1)
  # Each event has its own covered samples. RBPs with the same missingness share one projection.
  for c,s,ivs in clusters:
   keys=[(c,s,l,r) for l,r in ivs]
   if len(keys)<2: continue
   counts=np.array([[j.get(k,0) for k in keys] for j in js],dtype=float); totals=counts.sum(axis=1)
   for k,counts1 in zip(keys,counts.T):
    if k in annotated or np.count_nonzero(counts1[valid])<a.min_carriers: continue
    base=valid & (totals>=a.min_reads)
    if base.sum()<a.min_samples: continue
    patterns=defaultdict(list)
    for i in genes:
     mask=base & np.isfinite(logged[:,i]); patterns[np.packbits(mask).tobytes()].append(i)
    event_tested=False
    for group in patterns.values():
     mask=base & np.isfinite(logged[:,group[0]]); n=int(mask.sum()); xx=X[mask]
     rank=np.linalg.matrix_rank(xx); df=n-rank-1
     if n<a.min_samples or df<=1: continue
     y=counts1[mask]/totals[mask]; Q,_=np.linalg.qr(xx,mode='reduced') if rank==xx.shape[1] else (None,None)
     if Q is None:
      u,sv,_=np.linalg.svd(xx,full_matrices=False); Q=u[:,:rank]
     yr=y-Q@(Q.T@y); yy=float(yr@yr)
     if yy<1e-12: continue
     for start in range(0,len(group),a.batch_size):
      idx=group[start:start+a.batch_size]; raw=logged[np.ix_(mask,idx)]; xr=raw-Q@(Q.T@raw); ss=np.sum(xr*xr,axis=0); good=ss>1e-12
      if not np.any(good): continue
      ii=np.asarray(idx)[good]; corr=np.clip((xr[:,good].T@yr)/np.sqrt(ss[good]*yy),-1+1e-12,1-1e-12)
      pv=2*tdist.sf(np.abs(corr)*np.sqrt(df/(1-corr*corr)),df)
      jid=f'chr{c}:{s}:{k[2]}-{k[3]}'; carriers=int(np.count_nonzero(counts1[mask])); batch=[(geneids[int(i)],';'.join(sorted(gene_symbols[geneids[int(i)]])),jid,c,s,k[2],k[3],n,carriers,float(r),float(p),len(keys)) for i,r,p in zip(ii,corr,pv)]
      con.executemany('INSERT INTO pairs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,NULL)',batch); total_pairs+=len(batch); event_tested=True
    if event_tested: tested_events+=1
    if tested_events and tested_events%100==0: con.commit(); print('Events:',tested_events,'pairs:',total_pairs,flush=True)
  con.commit(); print('Calculating tissue-wide FDR for',total_pairs,'tests',flush=True)
  # The kth largest p-value has ascending rank m-k+1. Running minimum gives BH q-values, including ties.
  updates=[]; running=1.0
  for j,(rowid,pval) in enumerate(con.execute('SELECT rowid,p_value FROM pairs ORDER BY p_value DESC,rowid DESC')):
   running=min(running,float(pval)*total_pairs/(total_pairs-j)); updates.append((min(1.0,running),rowid))
   if len(updates)>=50000: con.executemany('UPDATE pairs SET fdr=? WHERE rowid=?',updates); updates=[]
  if updates: con.executemany('UPDATE pairs SET fdr=? WHERE rowid=?',updates)
  con.commit()
  with gzip.open(out/'all_tested_associations.tsv.gz','wt',newline='') as fh:
   writer=csv.writer(fh,delimiter='\t'); writer.writerow(cols+['fdr'])
   cur=con.execute('SELECT '+','.join(cols)+',fdr FROM pairs ORDER BY p_value')
   while True:
    rows=cur.fetchmany(50000)
    if not rows: break
    writer.writerows(rows)
  sig=pd.read_sql_query('SELECT '+','.join(cols)+',fdr FROM pairs WHERE fdr<=? ORDER BY fdr,p_value',con,params=(a.fdr,))
  sig.to_csv(out/'significant_edges.tsv',sep='\t',index=False)
  con.close(); dbpath.unlink()
  rbp=sig.groupby(['rbp_ensembl','rbp_gene'],as_index=False).agg(n_junctions=('junction_id','nunique'),junctions=('junction_id',lambda x:';'.join(sorted(set(x)))))
  ev=sig.groupby('junction_id',as_index=False).agg(n_rbps=('rbp_ensembl','nunique'),rbps=('rbp_gene',lambda x:';'.join(sorted(set(x)))))
  rbp.to_csv(out/'rbp_to_junctions.tsv',sep='\t',index=False); ev.to_csv(out/'junction_to_rbps.tsv',sep='\t',index=False)
  # Standalone interactive SVG: clicking a node highlights its connections.
  edges=sig.head(a.max_map_edges).to_dict('records'); data=json.dumps(edges).replace('</','<\\/')
  html='''<!doctype html><meta charset="utf-8"><title>RBP junction associations</title><style>body{font:14px sans-serif}svg{border:1px solid #ddd}line{stroke-opacity:.32}text{font-size:11px;cursor:pointer}input{margin:8px}</style><h2>RBP–unannotated junction associations</h2><p>Click a label to highlight its links. Red: positive; blue: negative. Display limited to the strongest edges; see tables for all results.</p><input id="q" placeholder="Search label"><svg id="g"></svg><script>const E=DATA,svg=document.querySelector('#g'),N='http://www.w3.org/2000/svg';let L=[...new Set(E.map(x=>x.rbp_gene))].sort(),R=[...new Set(E.map(x=>x.junction_id))].sort(),h=Math.max(400,24*Math.max(L.length,R.length)+60);svg.setAttribute('width',1200);svg.setAttribute('height',h);let lines=[];function el(n,a){let z=document.createElementNS(N,n);for(let k in a)z.setAttribute(k,a[k]);svg.appendChild(z);return z} E.forEach(e=>{let l=el('line',{x1:250,y1:30+24*L.indexOf(e.rbp_gene),x2:750,y2:30+24*R.indexOf(e.junction_id),stroke:e.partial_r>0?'#d45d57':'#447bb5','stroke-width':1.5});lines.push([l,e])}); function label(s,x,y){let t=el('text',{x,y});t.textContent=s;t.onclick=()=>lines.forEach(([l,e])=>l.style.strokeOpacity=(e.rbp_gene===s||e.junction_id===s)?1:.03)} L.forEach((s,i)=>label(s,10,34+24*i));R.forEach((s,i)=>label(s,760,34+24*i));document.querySelector('#q').oninput=e=>{let q=e.target.value.toLowerCase();lines.forEach(([l,x])=>l.style.strokeOpacity=!q||x.rbp_gene.toLowerCase().includes(q)||x.junction_id.toLowerCase().includes(q)? .8:.03)};</script>'''.replace('DATA',data)
  (out/'matching_map.html').write_text(html)
  (out/'analysis_summary.json').write_text(json.dumps(dict(cohort=cohort,tissue=tissue.name,samples=len(samples),complete_cases=int(valid.sum()),covariates=covused,tested_events=tested_events,tested_pairs=total_pairs,significant_pairs=len(sig),fdr_threshold=a.fdr,shared_cluster='junctions connected by shared splice site; pooled junction reads >= 5',coordinate_convention='junc end + 1; compare exon-end to next exon-start',map_edges=min(len(sig),a.max_map_edges)),indent=2))
  print(cohort,tissue.name,'tested pairs',total_pairs,'significant',len(sig),flush=True)
  marker.write_text(json.dumps({'settings':fingerprint,'tested_pairs':total_pairs,'significant_pairs':len(sig)}))
for cohort in (['target','answer'] if a.cohort=='both' else [a.cohort]): run(cohort)
print('Outputs:',a.outdir,flush=True)
PY
