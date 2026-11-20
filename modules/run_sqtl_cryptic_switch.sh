#!/bin/bash
#SBATCH --job-name=sqtl_cryptic_switch
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=sqtl_cryptic_switch.out
#SBATCH --error=sqtl_cryptic_switch.err

# Submit: sbatch run_sqtl_annotated_cryptic.sh --tissues Frontal_Cortex
# Omit --tissues to run all five Target ALS tissues. Additional options:
# --root PATH --outdir PATH --gtf PATH --min-n 5
set -euo pipefail
TASK_PYTHON=/home/zw529/donglab/pipelines/modules/miniconda3/bin/python
"$TASK_PYTHON" -I -u - "$@" <<'PYTHON'
#!/usr/bin/env python3
"""Rank opposite genotype-associated changes in connected annotated/cryptic junctions."""
import argparse, csv, re
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

TISSUES=['Cortex', 'Spinal_Cord', 'Cerebellum','Frontal_Cortex','Motor_Cortex','Cervical_Spinal_Cord','Lumbar_Spinal_Cord']

def jid(s):
    m=re.fullmatch(r'([^:]+):([+-]):(\d+)-(\d+)',str(s))
    if not m: raise ValueError(f'Invalid junction ID: {s}')
    return m.group(1),m.group(2),int(m.group(3)),int(m.group(4))

def scan_matrix(path, wanted):
    result={}
    with path.open() as f:
        header=f.readline().rstrip('\r\n').split('\t')[1:]
        if len(set(header)) != len(header): raise ValueError(f'Duplicate subject IDs: {path}')
        for line in f:
            parts=line.rstrip('\r\n').split('\t'); key=parts[0]
            if key not in wanted: continue
            if key in result: raise ValueError(f'Duplicate row {key}: {path}')
            if len(parts)!=len(header)+1: raise ValueError(f'Wrong row length: {key}')
            result[key]=pd.to_numeric(pd.Series(parts[1:]),errors='coerce').to_numpy(dtype=float)
    return header,result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.home()/'donglab/data/target_ALS')
    p.add_argument('--outdir',type=Path)
    p.add_argument('--gtf',type=Path,default=Path.home()/'donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf')
    p.add_argument('--tissues',nargs='+',default=TISSUES)
    p.add_argument('--min-n',type=int,default=5)
    a=p.parse_args(); out=a.outdir or a.root/'QTL/sqtl_annotated_cryptic'; out.mkdir(parents=True,exist_ok=True)
    if a.min_n < 1: p.error('--min-n must be positive')
    # Whole transcripts, sorted by genomic position; LeafViz uses exon end and next exon start.
    tx=defaultdict(list)
    with a.gtf.open() as f:
        for line in f:
            if line.startswith('#'): continue
            x=line.rstrip('\n').split('\t')
            if len(x)<9 or x[2]!='exon': continue
            m=re.search(r'transcript_id "([^"]+)"',x[8])
            if m: tx[(x[0],x[6],m.group(1))].append((int(x[3]),int(x[4])))
    annotated=set()
    for (c,strand,_),exons in tx.items():
        exons=sorted(set(exons))
        annotated.update((c,strand,left[1],right[0]) for left,right in zip(exons,exons[1:]) if left[1]<right[0])
    del tx
    print('Annotated junctions:',len(annotated),flush=True)
    audit=[]
    for tissue in a.tissues:
        q=a.root/tissue/'sQTL'; folder=out/tissue; folder.mkdir(parents=True,exist_ok=True)
        sig=pd.read_csv(q/'results'/f'{tissue}_sQTL.FDR0.05.txt',sep='\t')
        sig['FDR']=pd.to_numeric(sig['FDR'],errors='coerce')
        sig=sig[sig.FDR.between(0,0.05)].sort_values('FDR').drop_duplicates(['junction_id','snpid'])
        anchors=set(sig.junction_id); anchor_coords={j:jid(j) for j in anchors}
        starts=defaultdict(list); ends=defaultdict(list)
        for j,(c,s,b,e) in anchor_coords.items(): starts[(c,s,b)].append(j); ends[(c,s,e)].append(j)
        # Discover connected partners among the QTL matrix rows; do not infer zeros.
        partners=defaultdict(set)
        with (q/f'splicing_{tissue}.txt').open() as f:
            next(f)
            for line in f:
                j=line.split('\t',1)[0]; c,s,b,e=jid(j)
                for anchor in set(starts.get((c,s,b),[])+ends.get((c,s,e),[])):
                    if j!=anchor and ((c,s,b,e) in annotated) != (anchor_coords[anchor] in annotated): partners[anchor].add(j)
        wanted=anchors | {j for pp in partners.values() for j in pp}
        subjects,psis=scan_matrix(q/f'splicing_{tissue}.txt',wanted)
        gs,gts=scan_matrix(q/f'snp_{tissue}.txt',set(sig.snpid))
        if set(gs)!=set(subjects): raise ValueError(f'Subject sets differ: {tissue}')
        idx={s:i for i,s in enumerate(gs)}; order=[idx[s] for s in subjects]
        loc={}
        with (q/'snp_location.txt').open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                if r['snpid'] not in gts: continue
                k=(r['chr'].removeprefix('chr'),int(r['pos']))
                if r['snpid'] in loc and loc[r['snpid']]!=k: raise ValueError('Ambiguous SNP location '+r['snpid'])
                loc[r['snpid']]=k
        targets=set(loc.values()); alleles=defaultdict(set)
        with (a.root/'QTL/plink/joint_all_chrs_matrixEQTL.raw').open() as f:
            for token in f.readline().split()[6:]:
                m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([^:]+):([^_]+)_(.+)',token)
                if m:
                    c,pos,ref,alt,counted=m.groups(); k=(c,int(pos))
                    if k in targets: alleles[k].add((ref,alt,counted))
        records=[]; excluded=defaultdict(int)
        for row in sig.to_dict('records'):
            anchor=row['junction_id']; snp=row['snpid']; aa=alleles.get(loc.get(snp),set())
            if anchor not in psis or snp not in gts: excluded['missing_matrix_row']+=1; continue
            if len(aa)!=1: excluded['missing_or_ambiguous_alleles']+=1; continue
            ref,alt,counted=next(iter(aa))
            if counted not in (ref,alt) or ',' in alt or ref==alt: excluded['unsupported_alleles']+=1; continue
            g=gts[snp][order]
            if np.any(np.isfinite(g)&~np.isin(g,[0,1,2])): raise ValueError('Non-hardcall genotype: '+snp)
            g=2-g if counted==ref else g
            for partner in sorted(partners.get(anchor,())):
                an=anchor if anchor_coords[anchor] in annotated else partner
                cr=partner if an==anchor else anchor
                x,y=psis[an],psis[cr]
                if np.any(np.isfinite(x)&((x<0)|(x>1))) or np.any(np.isfinite(y)&((y<0)|(y>1))): raise ValueError('PSI outside [0,1]')
                valid=np.isfinite(x)&np.isfinite(y); groups=[valid&(g==i) for i in (0,1,2)]
                ns=[int(mask.sum()) for mask in groups]
                if ns[0]<a.min_n or ns[2]<a.min_n: excluded['insufficient_paired_homozygotes']+=1; continue
                xm=[100*float(x[mask].mean()) if n else np.nan for mask,n in zip(groups,ns)]
                ym=[100*float(y[mask].mean()) if n else np.nan for mask,n in zip(groups,ns)]
                dx=xm[2]-xm[0]; dy=ym[2]-ym[0]
                if dx*dy>=0: excluded['not_opposing']+=1; continue
                rec=dict(row); rec.update(tissue=tissue,annotated_junction=an,cryptic_junction=cr,partner_junction=partner,anchor_annotation='annotated' if an==anchor else 'cryptic',variant_chr='chr'+loc[snp][0],variant_pos=loc[snp][1],ref=ref,alt=alt,plink_counted_allele=counted,n_ref_ref=ns[0],n_het=ns[1],n_hom_alt=ns[2],annotated_ref_pct=xm[0],annotated_het_pct=xm[1],annotated_alt_pct=xm[2],cryptic_ref_pct=ym[0],cryptic_het_pct=ym[1],cryptic_alt_pct=ym[2],annotated_delta_pp=dx,cryptic_delta_pp=dy,switch_score_pp=min(abs(dx),abs(dy)))
                records.append(rec)
        df=pd.DataFrame(records)
        if not df.empty:
            df=df.sort_values(['switch_score_pp','FDR','junction_id','snpid','partner_junction'],ascending=[False,True,True,True,True]); df.insert(0,'rank',range(1,len(df)+1))
        else: df=pd.DataFrame(columns=['rank','tissue','junction_id','snpid','annotated_junction','cryptic_junction','switch_score_pp'])
        df.to_csv(folder/'ranked_switches.tsv',sep='\t',index=False,na_rep='NA')
        audit.append(dict(tissue=tissue,significant_pairs=len(sig),anchors_with_matrix_partners=len(partners),ranked_pairs=len(df),**excluded))
        print(tissue,'ranked',len(df),'exclusions',dict(excluded),flush=True)
    pd.DataFrame(audit).fillna(0).to_csv(out/'audit.tsv',sep='\t',index=False)
if __name__=='__main__': main()

PYTHON
