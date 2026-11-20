#!/usr/bin/env bash
#SBATCH --job-name=unc13a_stmn2
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/unc13a_stmn2_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/unc13a_stmn2_%j.err
set -euo pipefail

# Override these variables at submission if needed: sbatch --export=ALL,OUTDIR=/path ...
export TARGET_ROOT="${TARGET_ROOT:-$HOME/donglab/data/target_ALS}"
export ANSWER_ROOT="${ANSWER_ROOT:-$HOME/donglab/data/answer_ALS}"
# Results (tables and PDF) are written here; logs are written to QTL above.
export OUTDIR="${OUTDIR:-/home/zw529/donglab/data/target_ALS/QTL/UNC13A_STMN2}"
mkdir -p "$OUTDIR"
command -v Rscript >/dev/null || { echo 'ERROR: Rscript is required for plots; load your cluster R module before submitting.' >&2; exit 1; }
python3 - <<'PY'
import csv, os, sys, math, statistics
from pathlib import Path
from collections import Counter, defaultdict

out=Path(os.environ['OUTDIR']); out.mkdir(parents=True,exist_ok=True)
roots={'targetALS':Path(os.environ['TARGET_ROOT']), 'answerALS':Path(os.environ['ANSWER_ROOT'])}
remap={}
for x in ['Motor Cortex Lateral','Motor Cortex Medial','Lateral Motor Cortex','Medial Motor Cortex','Primary Motor Cortex L','Primary Motor Cortex M','Cortex_Motor_Unspecified','Cortex_Motor_BA4','BA4 Motor Cortex','Lateral_motor_cortex','Motor_Cortex_Lateral','Motor_Cortex_Medial','Lateral_Motor_Cortex','Medial_Motor_Cortex','Primary_Motor_Cortex_L','Primary_Motor_Cortex_M','Motor_Cortex']: remap[x]='Motor_Cortex'
for x in ['Spinal_Cord_Cervical','Cervical Spinal Cord','Cervical_spinal_cord','Spinal_cord_Cervical','Cervical_Spinal_Cord']: remap[x]='Cervical_Spinal_Cord'
for x in ['Lumbar Spinal Cord','Spinal_Cord_Lumbosacral','Lumbosacral_Spinal_Cord','Lumbar_spinal_cord','Lumbar_Spinal_Cord']: remap[x]='Lumbar_Spinal_Cord'
remap.update({'Frontal Cortex':'Frontal_Cortex','Frontal_Cortex':'Frontal_Cortex','Cerebellum':'Cerebellum'})
# LeafCutter .junc is BED-like: start is zero-based, end is exclusive. The observed PSI equivalent of .junc 17642591-17642844 is 17642591-17642845.
# Exact junctions, not +/- 1 matches. Shared boundary supports either isoform.
events={'shared':(17641556,17642413),'exon_128':(17642541,17642844),'exon_178':(17642591,17642844),'skipping':(17641556,17642844),'stmn2_cryptic':(79611214,79616821),'stmn2_canonical':(79611214,79636801)}
rows=[]; matches=[]; issues=[]
for dataset,root in roots.items():
    if not root.is_dir(): issues.append([dataset,str(root),'dataset root missing']); continue
    for processed in sorted(root.glob('*/RNAseq/Processed')):
        tissue=processed.parent.parent.name
        if dataset=='targetALS' and tissue not in remap: continue
        if dataset=='answerALS' and tissue!='iPSC_Motor_Neuron': continue
        mapped=remap.get(tissue,tissue)
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
                            if len(a)<6: continue
                            chrom=a[0][3:] if a[0].lower().startswith('chr') else a[0]
                            if (chrom,a[5]) not in [('19','-'),('8','+')]: continue
                            try: start,end,reads=int(a[1]),int(a[2]),int(a[4])
                            except ValueError: bad+=1; continue
                            # Preserve all nearby endpoint hits for convention auditing.
                            if (chrom=='19' and any(abs(x-y)<=1 for x in (start,end) for y in (17642414,17642541,17642591))) or (chrom=='8' and any(abs(x-y)<=1 for x in (start,end) for y in (79616822,79636802))):
                                matches.append([dataset,mapped,tissue,sd.name,a[0],a[5],start,end,reads,str(paths[0])])
                            for k,coord in events.items():
                                if (start,end)==coord and ((k.startswith('stmn2_') and chrom=='8') or (not k.startswith('stmn2_') and chrom=='19')): sums[k]+=reads
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
print('Output:',out,'; samples:',len(rows),'; issues:',len(issues),flush=True)
if not rows: sys.exit('No sample directories found; check dataset roots.')
PY

Rscript - "$OUTDIR" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
dir <- args[1]
d <- read.delim(file.path(dir, 'sample_junctions_STMN2.tsv'), check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('', 'NA'))
d <- d[grepl('^ok', d$junc_status), , drop=FALSE]
ratio <- function(a,b) ifelse(!is.na(a) & !is.na(b) & a+b>0, a/(a+b), NA_real_)
d$stmn2_usage <- ratio(d$stmn2_cryptic_reads, d$stmn2_canonical_reads)
d$unc13a_shared_usage <- ratio(d$shared_reads, d$skipping_reads)
den <- d$exon_128_reads + d$exon_178_reads + d$skipping_reads
d$unc13a_128_usage <- ifelse(!is.na(den) & den>0, d$exon_128_reads/den, NA_real_)
d$unc13a_178_usage <- ifelse(!is.na(den) & den>0, d$exon_178_reads/den, NA_real_)
write.table(d, file.path(dir,'sample_junction_usage.tsv'), sep='\t', row.names=FALSE, quote=FALSE, na='NA')
if (!nrow(d)) stop('No readable junction files: see issues.tsv')
pdf(file.path(dir,'UNC13A_STMN2_junction_usage.pdf'),width=12,height=4.5)
groups <- unique(d[c('dataset','tissue')]); stats <- list()
for (i in seq_len(nrow(groups))) {
 z <- d[d$dataset==groups$dataset[i] & d$tissue==groups$tissue[i],,drop=FALSE]
 par(mfrow=c(1,3),mar=c(4.5,4.5,3,1),oma=c(0,0,2,0))
 for (v in c('unc13a_shared_usage','unc13a_128_usage','unc13a_178_usage')) {
  ok <- is.finite(z$stmn2_usage) & is.finite(z[[v]]); x <- z$stmn2_usage[ok]; y <- z[[v]][ok]
  rho <- p <- NA_real_; status <- 'insufficient variation'
  if (length(x)>=3 && length(unique(x))>1 && length(unique(y))>1) {
   test <- suppressWarnings(cor.test(x,y,method='spearman',exact=FALSE)); rho <- unname(test$estimate); p <- test$p.value; status <- 'ok'
  }
  stats[[length(stats)+1]] <- data.frame(dataset=groups$dataset[i],tissue=groups$tissue[i],metric=v,n=length(x),rho=rho,p_value=p,status=status)
  plot(x,y,xlim=c(0,1),ylim=c(0,1),pch=16,col=adjustcolor('steelblue4',alpha.f=.45),xlab='STMN2 cryptic / (cryptic + canonical)',ylab='UNC13A junction usage',main=sub('unc13a_', 'UNC13A ', sub('_usage$', '', v)),cex.main=1.15)
  if (!is.na(rho)) { abline(lm(y ~ x), col='firebrick', lwd=1.5); legend('topright', legend=sprintf('Spearman rho = %.3f\nn = %d',rho,length(x)), bty='n', cex=.72) } else legend('topright', legend=paste('n =',length(x),'\nrho unavailable'), bty='n', cex=.72)
 }
 mtext(paste(groups$dataset[i],groups$tissue[i]),outer=TRUE,side=3,line=.5)
}
dev.off()
write.table(do.call(rbind,stats),file.path(dir,'spearman_usage_by_tissue.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
RSCRIPT
