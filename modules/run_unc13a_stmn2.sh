#!/usr/bin/env bash
#SBATCH --job-name=unc13a_stmn2
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/unc13a_stmn2.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/unc13a_stmn2.err
set -euo pipefail

# Override these variables at submission if needed: sbatch --export=ALL,OUTDIR=/path ...
export TARGET_ROOT="${TARGET_ROOT:-$HOME/donglab/data/target_ALS}"
export ANSWER_ROOT="${ANSWER_ROOT:-$HOME/donglab/data/answer_ALS}"
export OUTDIR="${OUTDIR:-$HOME/donglab/data/target_ALS/QTL/UNC13A_STMN2}"
mkdir -p "$OUTDIR"
python3 - <<'PY'
import csv, os, sys, math, statistics
from pathlib import Path
from collections import Counter, defaultdict

out=Path(os.environ['OUTDIR']); out.mkdir(parents=True,exist_ok=True)
roots={'targetALS':Path(os.environ['TARGET_ROOT']), 'answerALS':Path(os.environ['ANSWER_ROOT'])}
remap={}
for x in ['Motor Cortex Lateral','Motor Cortex Medial','Lateral Motor Cortex','Medial Motor Cortex','Primary Motor Cortex L','Primary Motor Cortex M','Cortex_Motor_Unspecified','Cortex_Motor_BA4','BA4 Motor Cortex','Lateral_motor_cortex','Motor_Cortex_Lateral','Motor_Cortex_Medial','Lateral_Motor_Cortex','Medial_Motor_Cortex','Primary_Motor_Cortex_L','Primary_Motor_Cortex_M']: remap[x]='Motor_Cortex'
for x in ['Spinal_Cord_Cervical','Cervical Spinal Cord','Cervical_spinal_cord','Spinal_cord_Cervical']: remap[x]='Cervical_Spinal_Cord'
for x in ['Lumbar Spinal Cord','Spinal_Cord_Lumbosacral','Lumbosacral_Spinal_Cord','Lumbar_spinal_cord']: remap[x]='Lumbar_Spinal_Cord'
remap.update({'Frontal Cortex':'Frontal_Cortex','Cerebellum':'Cerebellum'})
# LeafCutter .junc is BED-like: start is zero-based, end is exclusive. The observed PSI equivalent of .junc 17642591-17642844 is 17642591-17642845.
# Exact junctions, not +/- 1 matches. Shared boundary supports either isoform.
events={'shared':(17641556,17642413),'exon_128':(17642541,17642844),'exon_178':(17642591,17642844),'skipping':(17641556,17642844),'stmn2_cryptic':(79611214,79616821),'stmn2_canonical':(79611214,79636801)}
rows=[]; matches=[]; issues=[]
for dataset,root in roots.items():
    if not root.is_dir(): issues.append([dataset,str(root),'dataset root missing']); continue
    for processed in sorted(root.glob('*/RNAseq/Processed')):
        tissue=processed.parent.parent.name; mapped=remap.get(tissue,tissue)
        for sd in sorted(x for x in processed.iterdir() if x.is_dir()):
            base={'dataset':dataset,'tissue':mapped,'original_tissue':tissue,'sample':sd.name}
            paths=sorted((sd/'leafcutter'/'juncs').glob('*.junc'))
            row=dict(base,junc_file='',junc_status='')
            for key in events: row[key+'_reads']=''
            if len(paths)!=1:
                row['junc_status']='missing' if not paths else 'multiple_files'; issues.append([dataset,str(sd),row['junc_status']])
            else:
                row['junc_file']=str(paths[0]); sums=Counter(); bad=0
                try:
                    with paths[0].open() as f:
                        for line in f:
                            a=line.split()
                            if len(a)<6 or (a[0],a[5]) not in [('chr19','-'),('chr8','+')]: continue
                            try: start,end,reads=int(a[1]),int(a[2]),int(a[4])
                            except ValueError: bad+=1; continue
                            # Preserve all nearby endpoint hits for convention auditing.
                            if (a[0]=='chr19' and any(abs(x-y)<=1 for x in (start,end) for y in (17642414,17642541,17642591))) or (a[0]=='chr8' and any(abs(x-y)<=1 for x in (start,end) for y in (79616822,79636802))):
                                matches.append([dataset,mapped,tissue,sd.name,a[0],a[5],start,end,reads,str(paths[0])])
                            for k,coord in events.items():
                                if (start,end)==coord and ((k.startswith('stmn2_') and a[0]=='chr8') or (not k.startswith('stmn2_') and a[0]=='chr19')): sums[k]+=reads
                    row['junc_status']='ok' if not bad else 'ok_with_malformed_chr19_rows'
                    for k in events: row[k+'_reads']=sums[k]
                except (OSError,UnicodeError) as e:
                    row['junc_status']='read_error'; issues.append([dataset,str(paths[0]),str(e)])
            rows.append(row)
    print(dataset,': processed',sum(r['dataset']==dataset for r in rows),'sample directories',flush=True)

def write(name,fields,data):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(data)
fields=['dataset','tissue','original_tissue','sample','junc_status','shared_reads','exon_128_reads','exon_178_reads','skipping_reads','stmn2_cryptic_reads','stmn2_canonical_reads','junc_file']
write('sample_junctions_STMN2.tsv',fields,rows)
with (out/'nearby_boundary_audit.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['dataset','tissue','original_tissue','sample','chrom','strand','junc_start','junc_end','reads','file']); w.writerows(matches)
with (out/'issues.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['dataset','path','issue']); w.writerows(issues)
groups=defaultdict(list)
for r in rows: groups[(r['dataset'],r['tissue'])].append(r)
summary=[]
for (dataset,tissue),rr in sorted(groups.items()):
    ok=[r for r in rr if r['junc_status'].startswith('ok')]; paired=ok
    d=dict(dataset=dataset,tissue=tissue,samples=len(rr),junc_readable=len(ok),paired=len(paired))
    for k in events: d[k+'_positive']=sum(int(r[k+'_reads'])>0 for r in ok)
    summary.append(d)
write('tissue_summary.tsv',['dataset','tissue','samples','junc_readable','paired']+[k+'_positive' for k in events],summary)
# Compare junction reads within tissue, retaining zero counts from readable files.
try:
    from scipy.stats import spearmanr
except ImportError: spearmanr=None
stats=[]
for (dataset,tissue),rr in sorted(groups.items()):
    for k in ['shared','exon_128','exon_178']:
        pairs=[]
        for r in rr:
            if not r['junc_status'].startswith('ok'): continue
            try:
                x=float(r['stmn2_cryptic_reads']); y=float(r[k+'_reads'])
                if math.isfinite(x) and math.isfinite(y): pairs.append((x,y))
            except (ValueError,TypeError): pass
        d=dict(dataset=dataset,tissue=tissue,junction=k,n=len(pairs),positive=sum(y>0 for x,y in pairs),spearman_rho='',p_value='',status='')
        if spearmanr is None: d['status']='scipy_unavailable'
        elif len(pairs)<3 or len({x for x,y in pairs})<2 or len({y for x,y in pairs})<2: d['status']='insufficient_variation'
        else:
            z=spearmanr([x for x,y in pairs],[y for x,y in pairs]); d.update(spearman_rho=z.statistic,p_value=z.pvalue,status='ok')
        stats.append(d)
write('spearman_by_tissue.tsv',['dataset','tissue','junction','n','positive','spearman_rho','p_value','status'],stats)
try:
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for (dataset,tissue),rr in sorted(groups.items()):
        fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        for ax,k in zip(axes,['shared','exon_128','exon_178']):
            pairs=[]
            for r in rr:
                if not r['junc_status'].startswith('ok'): continue
                try:
                    x,y=float(r['stmn2_cryptic_reads']),float(r[k+'_reads'])
                    if math.isfinite(x) and math.isfinite(y): pairs.append((x,y))
                except (ValueError,TypeError): pass
            ax.scatter([x for x,y in pairs],[y for x,y in pairs],s=12,alpha=.45); ax.set_title(k.replace('_',' ')+' (n='+str(len(pairs))+')'); ax.set_xlabel('STMN2 cryptic junction reads'); ax.set_ylabel('Junction reads')
        fig.suptitle(dataset+' | '+tissue); fig.savefig(out/(dataset+'_'+tissue.replace('/','_')+'_STMN2.png'),dpi=160); plt.close(fig)
except ImportError: print('matplotlib unavailable; plots skipped',file=sys.stderr)
print('Output:',out,'; samples:',len(rows),'; issues:',len(issues),flush=True)
if not rows: sys.exit('No sample directories found; check dataset roots.')
PY
