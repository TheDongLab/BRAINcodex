#!/bin/bash
#SBATCH --job-name=sqtl_delta_psi
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/sqtl_splicing_difference/ranking_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/sqtl_splicing_difference/ranking_%j.err
# Usage: sbatch ~/donglab/pipelines/scripts/QTL/run_sQTL_splicing_difference.sh
# Check environment: bash run_sQTL_splicing_difference.sh --check
# Ranking: absolute difference of observed mean matrix PSI, ALT/ALT minus REF/REF.
# NA stays missing. Require >=5 usable subjects in BOTH homozygote groups.
# Plot up to 10 qualifying SNP-junction pairs per tissue, at most 2 SNPs per junction.
# Uses existing _sQTL_boxplot.R plus full-cluster PSI LeafViz plots.
# Ranking still uses QTL PSI; cluster arcs use unfiltered cluster counts.
# --tissues Cerebellum (uses existing per-sample refined cluster files automatically)
set -euo pipefail
TASK_PYTHON=/home/zw529/donglab/pipelines/modules/miniconda3/bin/python
TASK_SCRIPTS=/home/zw529/donglab/pipelines/scripts/QTL
TASK_RANK_ONLY=0
for task_arg in "$@"; do
    if [[ "$task_arg" == "--no-plots" ]]; then TASK_RANK_ONLY=1; fi
done
if [[ "$TASK_RANK_ONLY" == 0 ]]; then module load R; fi
# -I prevents mixed Python environments from loaded modules.
"$TASK_PYTHON" -I -c 'import sys,numpy,pandas; print("Python:",sys.executable); print("NumPy:",numpy.__file__); print("pandas:",pandas.__file__)'
if [[ "$TASK_RANK_ONLY" == 0 ]]; then
Rscript -e 'stopifnot(requireNamespace("data.table",quietly=TRUE),requireNamespace("ggplot2",quietly=TRUE),requireNamespace("gridExtra",quietly=TRUE),requireNamespace("svglite",quietly=TRUE)); cat("R plotting dependencies OK\n")'
for task_file in _sQTL_boxplot.R plot_sQTL_leafviz_PSI.R; do
    test -r "$TASK_SCRIPTS/$task_file"
done
fi
if [[ "${1:-}" == "--check" ]]; then exit 0; fi
"$TASK_PYTHON" -I -u - --scripts-dir "$TASK_SCRIPTS" "$@" <<'PYTHON'
#!/usr/bin/env python3
"""Rank significant sQTLs by observed homozygote mean PSI difference.

Stream large matrices; never impute missing PSI or assume PLINK counts ALT.
Ambiguous coordinate-to-allele mappings are excluded and reported.
"""
import argparse, csv, json, re, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd

TISSUES = ['Cerebellum','Frontal_Cortex','Motor_Cortex','Cervical_Spinal_Cord','Lumbar_Spinal_Cord']
PATTERNS = {'Cerebellum':'Cerebellum','Frontal_Cortex':'Frontal',
            'Motor_Cortex':'Motor.*Cortex|Cortex.*Motor|BA4',
            'Cervical_Spinal_Cord':'Cervical','Lumbar_Spinal_Cord':'Lumbar|Lumbosacral',
            'iPSC_Motor_Neuron':r'iPSC[ _-]*Motor[ _-]*Neuron'}

def matrix(path, wanted):
    data = {}
    with open(path) as f:
        samples = f.readline().split()[1:]
        if len(samples) != len(set(samples)): raise ValueError(f'Duplicate subjects: {path}')
        for line in f:
            key, _, tail = line.rstrip('\n').partition('\t')
            if key not in wanted: continue
            if key in data: raise ValueError(f'Duplicate row {key}: {path}')
            v = np.array([float(x) if x not in ('NA','NaN','nan','.','') else np.nan for x in tail.split('\t')])
            if len(v) != len(samples): raise ValueError(f'Wrong row length: {key}')
            data[key] = v
    return samples, data

def table(df, path):
    df.to_csv(path, sep='\t', index=False, na_rep='NA')

def rank_qualifying(rows):
    keep=(rows['n_ref_ref']>=5)&(rows['n_hom_alt']>=5)&rows['abs_delta_pp'].notna()
    ranked=rows.loc[keep].sort_values(['abs_delta_pp','FDR','junction_id','snpid'],ascending=[False,True,True,True]).reset_index(drop=True)
    ranked.insert(0,'rank',np.arange(1,len(ranked)+1))
    return ranked

def answer_metadata(meta,tissue,subjects,subject_column=None,sample_column=None):
    if meta.empty:raise ValueError('Answer ALS metadata is empty')
    m=meta.copy()
    subject_column=subject_column or ('externalsubjectid' if 'externalsubjectid' in m else m.columns[0])
    sample_column=sample_column or ('externalsampleid' if 'externalsampleid' in m else subject_column)
    for col in (subject_column,sample_column):
        if col not in m:raise ValueError(f'Missing metadata column: {col}')
    if 'tissue' in m:
        m=m[m.tissue.str.contains(PATTERNS[tissue],case=False,regex=True,na=False)]
    result=pd.DataFrame({'externalsubjectid':m[subject_column].astype(str).str.strip(),
                         'externalsampleid':m[sample_column].astype(str).str.strip()})
    result=result[result.externalsubjectid.isin(subjects)].drop_duplicates()
    if result.empty:raise ValueError('No Answer ALS metadata IDs match matrix subjects; supply --sample-map or ID-column options')
    if result.externalsampleid.eq('').any():raise ValueError('Missing Answer ALS RNA sample IDs')
    if result.externalsubjectid.duplicated().any() or result.externalsampleid.duplicated().any():
        raise ValueError('Ambiguous donor/sample mapping. Supply --sample-map with one RNA sample per matrix subject, matching the samples used to build the PSI matrix.')
    result['tissue']=tissue
    # The input matrices define the iPSC cohort; no postmortem/RIN filtering.
    return result

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--experiment',choices=['target_ALS','answer_ALS'],default=os.environ.get('EXPERIMENT','target_ALS'))
    ap.add_argument('--root',type=Path)
    ap.add_argument('--metadata',type=Path,help='Answer ALS defaults to aals_dataportal_datatable.csv')
    ap.add_argument('--sample-map',type=Path,help='Answer ALS TSV: externalsubjectid, externalsampleid; optional tissue')
    ap.add_argument('--subject-column',help='Answer ALS metadata column matching QTL matrix subjects')
    ap.add_argument('--sample-column',help='Answer ALS metadata column identifying RNA samples')
    ap.add_argument('--outdir',type=Path)
    ap.add_argument('--scripts-dir',type=Path,default=Path.home()/'donglab/pipelines/scripts/QTL')
    ap.add_argument('--tissues', nargs='+', choices=TISSUES+['iPSC_Motor_Neuron'], default=None)
    ap.add_argument('--cluster-counts', help='Unfiltered LeafCutter counts; use {tissue} for multiple tissues')
    ap.add_argument('--no-plots',action='store_true',help='Calculate rankings only')
    ap.add_argument('--gtf',type=Path,default=Path.home()/'donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf')
    a=ap.parse_args()
    if a.experiment not in ('target_ALS','answer_ALS'):ap.error('Invalid EXPERIMENT environment variable')
    a.root=(a.root or Path.home()/'donglab/data'/a.experiment).expanduser()
    a.outdir=(a.outdir or a.root/'QTL/sqtl_splicing_difference').expanduser()
    a.tissues=a.tissues or (['iPSC_Motor_Neuron'] if a.experiment=='answer_ALS' else TISSUES)
    if a.experiment=='target_ALS' and any([a.sample_map,a.subject_column,a.sample_column]):
        ap.error('Sample-map and ID-column options require --experiment answer_ALS')
    if a.sample_map and any([a.metadata,a.subject_column,a.sample_column]):
        ap.error('Use --sample-map alone or --metadata with optional ID columns')
    if a.experiment=='answer_ALS' and not a.no_plots:
        boxplot=a.scripts_dir/'_sQTL_boxplot.R'
        if 'SQTL_OPTIONAL_RAW_FILE' not in boxplot.read_text():
            ap.error('Install the accompanying updated _sQTL_boxplot.R first; it accepts the cohort PLINK file as argument 8.')
    # The existing boxplot reads EXPERIMENT, not an eighth positional argument.
    # Pass the resolved CLI/environment choice explicitly to every R child.
    plot_env=os.environ.copy()
    plot_env['EXPERIMENT']=a.experiment
    if not a.no_plots:
        expected_raw=Path.home()/'donglab/data'/a.experiment/'QTL/plink/joint_all_chrs_matrixEQTL.raw'
        ranking_raw=a.root/'QTL/plink/joint_all_chrs_matrixEQTL.raw'
        if ranking_raw.resolve()!=expected_raw.resolve():
            try: same_raw=ranking_raw.samefile(expected_raw)
            except OSError: same_raw=False
            if not same_raw:
                ap.error('Your boxplot resolves PLINK under ~/donglab/data/EXPERIMENT. '
                         'For plotting, --root must use that same PLINK file; otherwise use --no-plots.')
    a.outdir.mkdir(parents=True,exist_ok=True)
    print('Experiment:',a.experiment,'Root:',a.root,'Tissues:',a.tissues,flush=True)
    if not a.no_plots and a.cluster_counts:
        if len(a.tissues)>1 and '{tissue}' not in a.cluster_counts:
            ap.error('Use --tissues for one tissue or a --cluster-counts path containing {tissue}')
        for t in a.tissues:
            if not Path(a.cluster_counts.replace('{tissue}',t)).is_file():
                ap.error(f'Cluster counts missing for {t}')
    sigs={}; locs={}; targets=set()
    for t in a.tissues:
        q=a.root/t/'sQTL'
        s=pd.read_csv(q/'results'/f'{t}_sQTL.FDR0.05.txt',sep='\t')
        s['FDR']=pd.to_numeric(s.FDR,errors='coerce')
        s=s[s.FDR.between(0,0.05)].sort_values('FDR').drop_duplicates(['junction_id','snpid'])
        sigs[t]=s; wanted=set(s.snpid); loc={}
        with (q/'snp_location.txt').open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['snpid'] in wanted:
                    k=(r['chr'].removeprefix('chr'),int(r['pos']))
                    if r['snpid'] in loc and loc[r['snpid']]!=k: raise ValueError('Ambiguous SNP location')
                    loc[r['snpid']]=k; targets.add(k)
        locs[t]=loc
        print(t,'significant pairs',len(s),'SNPs',len(wanted),flush=True)
    alleles={}
    raw=a.root/'QTL/plink/joint_all_chrs_matrixEQTL.raw'
    with raw.open() as f:
        for token in f.readline().split()[6:]:
            m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([^:]+):([^_]+)_(.+)',token)
            if not m: continue
            c,p,ref,alt,counted=m.groups(); k=(c,int(p))
            if k in targets: alleles.setdefault(k,set()).add((ref,alt,counted))
    print('Allele coordinates mapped',len(alleles),flush=True)
    meta=None
    if not a.no_plots:
        metadata_path=(a.sample_map or a.metadata or a.root/('aals_dataportal_datatable.csv' if a.experiment=='answer_ALS' else 'targetALS_rnaseq_metadata.csv')).expanduser()
        if a.experiment=='answer_ALS':
            meta=pd.read_csv(metadata_path,sep='\t' if metadata_path.suffix.lower() in ('.txt','.tsv') else ',',dtype=str,keep_default_na=False,encoding='utf-8-sig')
            meta.columns=meta.columns.str.strip()
            if a.sample_map and not {'externalsubjectid','externalsampleid'}.issubset(meta.columns):
                ap.error('--sample-map requires externalsubjectid and externalsampleid')
        else:meta=pd.read_csv(metadata_path,low_memory=False)
        print('Metadata:',metadata_path,flush=True)
    genes={}
    with a.gtf.open() as f:
        for line in f:
            if line.startswith('#'): continue
            z=line.split('\t')
            if len(z)<9 or z[2]!='gene': continue
            m=re.search(r'gene_name "([^"]+)"',z[8])
            if m: genes.setdefault((z[0],z[6]),[]).append((int(z[3]),int(z[4]),m[1]))
    tops=[]; audit=[]; plot_jobs=[]
    for t in a.tissues:
        print('Loading selected matrix rows:',t,flush=True)
        q=a.root/t/'sQTL'; out=a.outdir/t; out.mkdir(exist_ok=True)
        s=sigs[t]; subjects,psi=matrix(q/f'splicing_{t}.txt',set(s.junction_id))
        gs,gts=matrix(q/f'snp_{t}.txt',set(s.snpid))
        if set(gs)!=set(subjects): raise ValueError('PSI/genotype subject sets differ')
        gi={x:i for i,x in enumerate(gs)}; order=[gi[x] for x in subjects]
        records=[]; excluded=[]; mapped={}
        for row in s.to_dict('records'):
            snp=row['snpid']; j=row['junction_id']; k=locs[t].get(snp)
            aa=alleles.get(k,set())
            reason=None
            if j not in psi or snp not in gts: reason='missing_matrix_row'
            elif len(aa)!=1: reason='missing_or_ambiguous_alleles'
            else:
                ref,alt,counted=next(iter(aa))
                if counted not in (ref,alt) or ',' in alt or ref==alt: reason='unsupported_alleles'
            if reason:
                excluded.append(dict(row,reason=reason));continue
            x=psi[j]; g=gts[snp][order]
            if np.any(np.isfinite(x)&((x<0)|(x>1))): raise ValueError(f'Non-PSI value: {j}')
            if np.any(np.isfinite(g)&~np.isin(g,[0,1,2])): raise ValueError(f'Non-hardcall genotype: {snp}')
            altg=2-g if counted==ref else g
            mapped[snp]=(altg,ref,alt,counted,k)
            r=dict(row,tissue=t,variant_chr='chr'+k[0],variant_pos=k[1],ref=ref,alt=alt,plink_counted_allele=counted)
            for value,label in [(0,'ref_ref'),(1,'het'),(2,'hom_alt')]:
                v=x[(altg==value)&np.isfinite(x)]
                r['n_'+label]=len(v);r['mean_'+label+'_pct']=100*float(v.mean()) if len(v) else np.nan
            r['delta_alt_minus_ref_pp']=r['mean_hom_alt_pct']-r['mean_ref_ref_pct']
            r['abs_delta_pp']=abs(r['delta_alt_minus_ref_pp'])
            r['both_homozygotes_n_ge_5']=r['n_ref_ref']>=5 and r['n_hom_alt']>=5
            c,strand,span=j.split(':'); start,end=map(int,span.split('-'))
            r['gene_name']=';'.join(sorted({name for lo,hi,name in genes.get((c,strand),[]) if lo<=end and hi>=start})) or 'unannotated'
            records.append(r)
        allrows=pd.DataFrame(records)
        if allrows.empty: raise ValueError(f'No mapped pairs in {t}')
        rank=rank_qualifying(allrows)
        selected_pairs=rank.loc[rank.groupby('junction_id',sort=False).cumcount()<2].head(10)
        table(rank,out/'ranked_qualifying.tsv')
        table(rank.drop_duplicates('junction_id'),out/'best_snp_per_junction.tsv')
        rejected=allrows[~allrows.both_homozygotes_n_ge_5].copy()
        rejected['reason']='fewer_than_5_usable_subjects_in_one_or_both_homozygote_groups'
        table(rejected,out/'excluded_homozygote_counts.tsv')
        table(pd.DataFrame(excluded,columns=list(s.columns)+['reason']),out/'excluded_mapping.tsv')
        audit.append(dict(tissue=t,significant_pairs=len(s),qualifying_pairs=len(rank),selected_for_plots=len(selected_pairs),excluded_homozygote_counts=len(rejected),mapping_excluded=len(excluded)))
        mt=None
        if not a.no_plots:
            if a.experiment=='answer_ALS':
                mt=answer_metadata(meta,t,subjects,a.subject_column,a.sample_column)
            else:
                mt=meta.copy(); mismatch=a.root/t/'eQTL/potential_sex_mismatches.txt'
                if mismatch.exists(): mt=mt[~mt.externalsampleid.isin(mismatch.read_text().splitlines())]
                mt=mt[mt.externalsubjectid.isin(subjects)]
                mt=mt[mt.apply(lambda r:r.astype(str).str.contains(PATTERNS[t],case=False,regex=True).any(),axis=1)]
                mt['RIN_score']=mt.iloc[:,[16,17]].apply(pd.to_numeric,errors='coerce').fillna(0).max(axis=1)
                mt=mt[(pd.to_numeric(mt.post_mortem_interval_in_hours,errors='coerce')<=40)&(mt.RIN_score>=3)]
                mt=mt.sort_values('RIN_score',ascending=False).drop_duplicates('externalsubjectid')
        picks={}
        for top in selected_pairs.to_dict('records'):
            assert top['n_ref_ref']>=5 and top['n_hom_alt']>=5
            tops.append(top)
            picks[(top['junction_id'],top['snpid'])]=top
        # Small aligned exports let the existing boxplot script avoid reloading GB-sized matrices.
        pair_df=pd.DataFrame(list(picks),columns=['junction_id','snpid'])
        pair_df.to_csv(out/'top_pairs.tsv',sep='\t',index=False,header=False)
        for name,ids,values,label in [('top_splicing_matrix.tsv',list(dict.fromkeys(j for j,snp in picks)),psi,'geneid'),
                                      ('top_snp_matrix.tsv',list(dict.fromkeys(snp for j,snp in picks)),gts,'snpid')]:
            rows=[values[key] if label=='geneid' else values[key][order] for key in ids]
            frame=pd.DataFrame(rows,index=ids,columns=subjects);frame.index.name=label
            frame.to_csv(out/name,sep='\t',na_rep='NA')
        if a.no_plots:
            print(t,'qualifying',len(rank),'selected',len(selected_pairs),flush=True)
            continue
        cov=pd.read_csv(q/f'covariates_{t}_encoded.txt',sep='\t',index_col=0)
        cov[subjects].to_csv(out/'top_covariates.tsv',sep='\t',na_rep='NA')
        with (out/'top_snp_locations.tsv').open('w') as f:
            writer=csv.writer(f,delimiter='\t');writer.writerow(['snpid','chr','pos'])
            for snp in dict.fromkeys(snp for j,snp in picks):
                c,pos=locs[t][snp];writer.writerow([snp,'chr'+c,pos])
        for (j,snp),top in picks.items():
            prefix=re.sub(r'[^A-Za-z0-9_.+-]','_',t+'__'+j+'__'+snp)
            x=psi[j]; altg,ref,alt,counted,k=mapped[snp]
            ok=np.isfinite(x)&np.isfinite(altg)
            sample=pd.DataFrame({'externalsubjectid':subjects,'PSI':x,'ALT_dosage':altg})[ok]
            sample['GT']=sample.ALT_dosage.map({0:'0/0',1:'0/1',2:'1/1'})
            table(sample,out/(prefix+'.matrix_subjects.tsv'))
            table(sample[['externalsubjectid','GT']],out/(prefix+'.genotypes.tsv'))
            selected=mt[mt.externalsubjectid.isin(sample.externalsubjectid)]
            selected.to_csv(out/(prefix+'.metadata.csv'),index=False)
            if set(selected.externalsubjectid)!=set(sample.externalsubjectid): raise ValueError(f'Sample metadata mismatch for {prefix}')
            plot_jobs.append(dict(top,prefix=prefix,outdir=str(out),metadata=str(out/(prefix+'.metadata.csv')),genotypes=str(out/(prefix+'.genotypes.tsv'))))
        print(t,'qualifying (both homozygotes n>=5)',len(rank),'top 10 (max 2 SNPs per junction)',selected_pairs[['junction_id','snpid','abs_delta_pp','n_ref_ref','n_hom_alt']].to_dict('records'),flush=True)
    table(pd.DataFrame(tops),a.outdir/'top10_per_tissue.tsv');table(pd.DataFrame(audit),a.outdir/'audit.tsv')
    (a.outdir/'plot_jobs.json').write_text(json.dumps(plot_jobs,indent=2))
    if not a.no_plots:
        for t in a.tissues:
            out=a.outdir/t
            if not (out/'top_pairs.tsv').stat().st_size: continue
                        cmd = [
                'Rscript',
                str(a.scripts_dir / '_sQTL_boxplot.R'),
                str(out / 'top_pairs.tsv'),
                str(out / 'top_snp_matrix.tsv'),
                str(out / 'top_splicing_matrix.tsv'),
                str(out / 'top_covariates.tsv'),
                str(out / 'top_snp_locations.tsv'),
                str(out),
                t,
            ]
            with (out / 'boxplots.log').open('w') as log:
                subprocess.run(
                    cmd,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    env=plot_env,
                )
        for p in plot_jobs:
            command=['Rscript',str(a.scripts_dir/'plot_sQTL_leafviz_PSI.R'),'--tissue',p['tissue'],'--tissue-regex',PATTERNS[p['tissue']],
                '--anchor',p['junction_id'],'--gene',p['gene_name'].split(';')[0],
                '--variant-chr',p['variant_chr'],'--variant-pos',str(p['variant_pos']),
                '--variant-id',p['snpid'],'--variant-ref',p['ref'],'--variant-alt',p['alt'],
                '--genotypes',p['genotypes'],'--metadata',p['metadata'],'--data-root',str(a.root),
                '--gtf',str(a.gtf),'--outdir',p['outdir'],'--prefix',p['prefix']]
            print('Plotting',p['prefix'],flush=True)
            for view in ['PSI']:
                paired=command.copy()
                paired[1]=str(a.scripts_dir/f'plot_sQTL_leafviz_{view}.R')
                if a.cluster_counts:
                    paired += ['--cluster-counts',a.cluster_counts.replace('{tissue}',p['tissue'])]
                with open(Path(p['outdir'])/(p['prefix']+f'.{view}.plot.log'),'w') as log:
                    subprocess.run(paired,stdout=log,stderr=subprocess.STDOUT,check=True)
                    subprocess.run(paired,stdout=log,stderr=subprocess.STDOUT,check=True,env=plot_env)
    print('COMPLETE',flush=True)

if __name__=='__main__': main()

PYTHON
