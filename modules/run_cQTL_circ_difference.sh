#!/bin/bash
#SBATCH --job-name=cqtl_delta_circ
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --time=1-23:00:00
#SBATCH --partition=week
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/cqtl_circ_difference_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/cqtl_circ_difference_%j.err

# Answer ALS: pass --experiment answer_ALS; override sbatch --output/--error for Answer ALS logs
# Optional --sample-map TSV: externalsubjectid, externalsampleid, optional tissue
# Usage for specific tissues: QTL_TISSUES="Cortex Spinal_Cord" sbatch --export=ALL run_cQTL_circ_difference.sh

set -euo pipefail
module load deepTools

command -v samtools >/dev/null
python -c 'import numpy,pandas,matplotlib,pyBigWig'
if [[ "${1:-}" == --check ]]; then exit 0; fi
python -u - "$@" <<'PYTHON'

# Strand-aware cQTL boundary figure; embedded to keep this job standalone.
import base64 as _cq64, gzip as _cqgz, types as _cqtypes
cq_boundary = _cqtypes.ModuleType('cqtl_boundary_plot')
exec(compile(_cqgz.decompress(_cq64.b64decode('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')), '<embedded cqtl_boundary_plot>', 'exec'), cq_boundary.__dict__)

import argparse, csv, json, re, subprocess, math, os, hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyBigWig
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FuncFormatter
CANONICAL_TISSUES=['Cerebellum','Frontal_Cortex','Motor_Cortex','Cervical_Spinal_Cord','Lumbar_Spinal_Cord','Cortex','Spinal_Cord','iPSC_Motor_Neuron']
GENOTYPE_LABELS={0:'Ref/Ref',1:'Het',2:'Hom Alt'}
TRACK_COLORS={0:'#228B22',1:'#4682B4',2:'#C75B7A'}
processed_root_cache={};processed_sample_cache={};primary_read_cache={};leafcutter_cache={}
def norm(s):return re.sub(r"[^a-z0-9]","",str(s).lower())

def chr_key(s):return re.sub(r"^chr","",str(s),flags=re.I).upper()

def safe_name(s):return re.sub(r"[^A-Za-z0-9._+-]+","_",str(s))[:180]

def tissue_display(x):return str(x).replace("_"," ")

def missing(v):return v is None or str(v).strip().lower() in {"","na","nan","none",".","null"}

def as_float(x):
    try:
        v=float(x);return v if np.isfinite(v) else None
    except:return None

def format_max_sig(v,sig=5):
    if v==0:return "0"
    d=max(0,sig-int(math.floor(math.log10(abs(v))))-1)
    return f"{v:.{d}f}".rstrip("0").rstrip(".")

def genomic_formatter(x,pos=None):return f"{format_max_sig(x/1e6)} Mb" if abs(x)>=1e6 else f"{format_max_sig(x/1e3)} kb" if abs(x)>=1e3 else format_max_sig(x)

def density_formatter(y,pos=None):return "0" if abs(y)<1e-12 else format_max_sig(abs(y),4)

def canonical_tissue(x):
    x=str(x).strip();return TISSUE_REMAP_NORM.get(norm(x),x)

def tissue_equal(a,b):
    if b=="Cortex":return bool(re.search(r"(?:^|[\W_])cortex(?:$|[\W_])",str(a),re.I))
    if b=="Spinal_Cord":return bool(re.search(r"(?:^|[\W_])spinal[\W_]*cord(?:$|[\W_])",str(a),re.I))
    return canonical_tissue(a)==canonical_tissue(b)

def parse_event_coordinates(event):
    """
    Parse the coordinate conventions used by the QTL outputs.

    cQTL example:
      chr10:100164003-100167406:-

    sQTL example:
      chr10:+:100381448-100394498

    Also accepts a few legacy colon-separated variants.
    """
    s=str(event).strip()

    patterns=[
        # sQTL: chr10:+:100381448-100394498
        (r"^(chr[^:]+):([+-]):(\d+)-(\d+)$", "strand_first"),

        # cQTL: chr10:100164003-100167406:-
        (r"^(chr[^:]+):(\d+)-(\d+):([+-])$", "strand_last"),

        # Legacy: chr10:100164003:100167406:-
        (r"^(chr[^:]+):(\d+):(\d+):([+-])$", "strand_last"),

        # No strand: chr10:100164003-100167406
        (r"^(chr[^:]+):(\d+)-(\d+)$", "no_strand"),

        # No strand: chr10:100164003:100167406
        (r"^(chr[^:]+):(\d+):(\d+)$", "no_strand")
    ]

    for pat,mode in patterns:
        m=re.match(pat,s)
        if not m:
            continue

        if mode=="strand_first":
            chrom,strand,start,end=m.groups()
        elif mode=="strand_last":
            chrom,start,end,strand=m.groups()
        else:
            chrom,start,end=m.groups()
            strand=None

        return {
            "chr":chrom,
            "start":int(start),
            "end":int(end),
            "strand":strand
        }

    return None

def processed_roots_for_tissue(ROOT,tissue):
    t=canonical_tissue(tissue)
    key=(str(ROOT),t,str(PROCESSED_ROOT))
    if key in processed_root_cache:return processed_root_cache[key]
    roots=[]
    if PROCESSED_ROOT is not None:
        roots=[PROCESSED_ROOT]
    else:
        preferred=ROOT/t/"RNAseq"/"Processed"
        if preferred.is_dir():roots.append(preferred)
        for p in ROOT.iterdir():
            if p.is_dir() and tissue_equal(p.name,t):
                q=p/"RNAseq"/"Processed"
                if q.is_dir() and q not in roots:roots.append(q)
    processed_root_cache[key]=roots
    print(f"  RNAseq roots for {t}: "+(" | ".join(map(str,roots)) if roots else "NONE"),flush=True)
    return roots

def resolve_sample_dir(ROOT,tissue,sid):
    key=(str(ROOT),canonical_tissue(tissue),sid,str(PROCESSED_ROOT))
    if key in processed_sample_cache:return processed_sample_cache[key]
    roots=processed_roots_for_tissue(ROOT,tissue)
    matches=set()
    for root in roots:
        for x in [sid,sid.replace("-","_"),sid.replace("_","-")]:
            p=root/x
            if p.is_dir():matches.add(p.resolve())
    # Keep the legacy broader matching for Target ALS only.
    if not matches and EXPERIMENT=='target_ALS':
        for root in roots:
            matches.update(p.resolve() for p in root.iterdir() if p.is_dir() and norm(p.name)==norm(sid))
    if len(matches)>1:
        raise ValueError(f'Ambiguous sample directory for {sid}: {sorted(map(str,matches))}')
    result=next(iter(matches),None)
    processed_sample_cache[key]=result
    return result

EXPERIMENT='target_ALS'
PROCESSED_ROOT=None

def primary_read_count(sd):
    bam=sd/"STAR.Aligned.sortedByCoord.out.bam";k=str(bam)
    if k in primary_read_cache:return primary_read_cache[k]
    if not bam.exists():primary_read_cache[k]=None;return None
    try:n=int(subprocess.check_output(["samtools","view","-F","0x100","-c",str(bam)],text=True).strip());primary_read_cache[k]=n if n>0 else None
    except:primary_read_cache[k]=None
    return primary_read_cache[k]

def rpm_value(reads,n):return float(reads)*1e6/n if n and n>0 else None

def bw_chrom(bw,chrom):
    cs=bw.chroms();return next((c for c in [chrom,chrom[3:] if chrom.startswith("chr") else "chr"+chrom] if c in cs),None)

def bigwig_strands(sd,chrom,start,end,n_bins):
    paths=[sd/"STAR.Aligned.sortedByCoord.out.plus.normalized.bw",sd/"STAR.Aligned.sortedByCoord.out.minus.normalized.bw"]
    def read(path):
        if not path.exists():return None
        bw=pyBigWig.open(str(path))
        try:
            c=bw_chrom(bw,chrom)
            if c is None:return None
            return np.abs(np.asarray([0 if v is None else float(v) for v in bw.stats(c,int(start),int(end),nBins=n_bins,type="mean")]))
        finally:bw.close()
    return read(paths[0]),read(paths[1]),str(paths[0]),str(paths[1])

def leafcutter_file(sd):return sd/"leafcutter"/"psi"/f"{sd.name}.leafcutter.PSI.tsv"

def load_leafcutter(path):
    k=str(path)
    if k in leafcutter_cache:return leafcutter_cache[k]
    if not path.exists():leafcutter_cache[k]=None;return None
    out={}
    try:
        with open(path) as f:
            for r in csv.DictReader(f,delimiter="\t"):
                try:key=(chr_key(r["chrom"]),r["strand"].strip(),int(float(r["start"])),int(float(r["end"])));out[key]=out.get(key,0)+float(r["reads"])
                except:continue
    except:out=None
    leafcutter_cache[k]=out;return out

def mean_junction_rpm(group,chrom,start,end,MIN_MEAN_JUNCTION_RPM,MAX_JUNCTIONS):
    sums=defaultdict(float);n=0
    for r in group:
        j=load_leafcutter(leafcutter_file(r["sample_dir"]));primary=r["primary_reads"]
        if j is None or not primary:continue
        n+=1
        for (c,st,s,e),reads in j.items():
            if c==chr_key(chrom) and s>=start and e<=end:sums[(st,s,e)]+=rpm_value(reads,primary)
    if not n:return {},0
    means={k:v/n for k,v in sums.items()};means={k:v for k,v in means.items() if v>=MIN_MEAN_JUNCTION_RPM}
    if MAX_JUNCTIONS>0 and len(means)>MAX_JUNCTIONS:
        keep=sorted(means,key=means.get,reverse=True)[:MAX_JUNCTIONS];means={k:means[k] for k in keep}
    return means,n

def parse_gtf_attrs(s):
    out={}
    for x in s.strip().split(";"):
        p=x.strip().split(" ",1)
        if len(p)==2:out[p[0]]=p[1].strip().strip('"')
    return out

def load_gene_index(gtf):
    # Match the zero-based, half-open convention used by the coverage plots.
    # GTF coordinates are one-based inclusive, so subtract one from the start.
    genes=defaultdict(list)
    with open(gtf) as f:
        for line in f:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)<9 or fields[2]!='gene':continue
            attrs=parse_gtf_attrs(fields[8])
            gene_id=attrs.get('gene_id','')
            gene_name=attrs.get('gene_name','')
            if not gene_id and not gene_name:continue
            genes[chr_key(fields[0])].append((int(fields[3])-1,int(fields[4]),fields[6],gene_name,gene_id))
    if not genes:raise ValueError(f'No gene annotations found in GTF: {gtf}')
    for chrom in genes:genes[chrom].sort(key=lambda row:row[0])
    return genes

def annotate_circ_gene(event,gene_index,cache):
    # Report all overlapping genes on the circRNA strand, not the SNP genes.
    # Multiple symbols/IDs are separated by semicolons; missing matches stay NA.
    if event in cache:return cache[event]
    coords=parse_event_coordinates(event)
    symbols=set();gene_ids=set()
    if coords is not None:
        start,end=sorted((coords['start'],coords['end']))
        strand=coords['strand']
        for left,right,gene_strand,symbol,gene_id in gene_index.get(chr_key(coords['chr']),[]):
            if left>=end:break
            if right<=start:continue
            if strand in ('+','-') and gene_strand!=strand:continue
            if symbol:symbols.add(symbol)
            if gene_id:gene_ids.add(gene_id)
    cache[event]={
        'gene_name':';'.join(sorted(symbols)) if symbols else None,
        'gene_id':';'.join(sorted(gene_ids)) if gene_ids else None
    }
    return cache[event]

def load_reference_exons(gtf,chrom,start,end):
    out=[]
    with open(gtf) as f:
        for l in f:
            if l.startswith("#"):continue
            x=l.rstrip().split("\t")
            if len(x)<9 or x[2]!="exon" or chr_key(x[0])!=chr_key(chrom):continue
            s,e=int(x[3])-1,int(x[4])
            if e<start or s>end:continue
            a=parse_gtf_attrs(x[8]);tid=a.get("transcript_id","")
            if not tid:continue
            out.append({"start":s,"end":e,"strand":x[6],"gene":a.get("gene_name",a.get("gene_id","Unknown")),"gene_id":a.get("gene_id",""),"transcript_id":tid,"transcript_name":a.get("transcript_name","")})
    return out

def draw_reference_track(ax,exons,start,end):
    transcripts=defaultdict(list);info={}
    for e in exons:
        tid=e["transcript_id"];transcripts[tid].append((e["start"],e["end"]))
        info[tid]={"gene":e["gene"],"strand":e["strand"],"transcript_name":e["transcript_name"]}
    if not transcripts:
        ax.text(.5,.45,"No GENCODE v49 transcripts",transform=ax.transAxes,ha="center",fontsize=9);ax.set_ylim(-.45,.55)
    else:
        order=sorted(transcripts,key=lambda tid:(min(s for s,e in transcripts[tid]),info[tid]["gene"],tid))
        for y,tid in enumerate(order):
            ints=sorted(transcripts[tid]);tx_start=max(start,min(s for s,e in ints));tx_end=min(end,max(e for s,e in ints))
            ax.plot([tx_start,tx_end],[y,y],color="black",lw=.7)
            for s,e in ints:
                s=max(s,start);e=min(e,end)
                if e>s:ax.add_patch(Rectangle((s,y-.14),e-s,.28,facecolor="black",edgecolor="black",lw=.35))
            gene=info[tid]["gene"];strand=info[tid]["strand"];tname=info[tid]["transcript_name"]
            label=f"{gene} | {tname} | {tid} ({strand})" if tname and tname!=tid else f"{gene} | {tid} ({strand})"
            ax.text(start,y+.19,label,fontsize=6.3,ha="left",va="bottom")
        ax.set_ylim(-.45,len(order)-.05)
    ax.set_xlim(start,end);ax.set_yticks([]);ax.set_ylabel("GENCODE v49",fontsize=12,labelpad=13)
    ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter));ax.locator_params(axis="x",nbins=5);ax.tick_params(axis="x",labelsize=11)
    for s in ["top","right","left"]:ax.spines[s].set_visible(False)

def draw_linear_arcs(ax,junctions,x,plus,minus,cmax,color):
    items=sorted(junctions.items(),key=lambda z:(z[0][0],(z[0][1]+z[0][2])/2,z[0][2]-z[0][1]));span=max(1,x[-1]-x[0])
    for i,((st,s,e),rpm) in enumerate(items):
        frac=min(1,max(0,(e-s)/span));mid=(s+e)/2;stagger=(i%3)*.10*cmax
        if st=="-":ly,ry=-float(np.interp(s,x,minus)),-float(np.interp(e,x,minus));base=min(ly,ry);apex=-cmax*(1.05+.42*math.sqrt(frac))-stagger
        else:ly,ry=float(np.interp(s,x,plus)),float(np.interp(e,x,plus));base=max(ly,ry);apex=cmax*(1.05+.42*math.sqrt(frac))+stagger
        path=MplPath([(s,ly),(mid,apex),(e,ry)],[MplPath.MOVETO,MplPath.CURVE3,MplPath.CURVE3]);lw=min(4,max(.7,math.log(rpm+1,2)))
        ax.add_patch(PathPatch(path,facecolor="none",edgecolor=color,lw=lw,alpha=.9))
        ax.text(mid,base+.58*(apex-base),f"{rpm:.3f}",ha="center",va="center",fontsize=7,bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=.45))

def draw_circ_arc(ax,start,end,rpm,cmax,color,strand):
    if rpm<=0:return
    outer=(-1 if strand=="-" else 1)*cmax*1.72
    path=MplPath([(start,0),(start,outer),(end,outer),(end,0)],[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
    lw=min(4,max(.7,math.log(rpm+1,2)));ax.add_patch(PathPatch(path,facecolor="none",edgecolor=color,lw=max(1.4,lw),linestyle="--",alpha=.95))
    ax.text((start+end)/2,outer*.63,f"{rpm:.3f}",ha="center",va="center",fontsize=7,bbox=dict(facecolor="white",edgecolor="none",alpha=.88,pad=.45))

TISSUE_REMAP={
"Motor Cortex Lateral":"Motor_Cortex","Motor Cortex Medial":"Motor_Cortex","Lateral Motor Cortex":"Motor_Cortex","Medial Motor Cortex":"Motor_Cortex",
"Primary Motor Cortex L":"Motor_Cortex","Primary Motor Cortex M":"Motor_Cortex","Cortex_Motor_Unspecified":"Motor_Cortex","Cortex_Motor_BA4":"Motor_Cortex",
"BA4 Motor Cortex":"Motor_Cortex","Lateral_motor_cortex":"Motor_Cortex","Motor_Cortex_Lateral":"Motor_Cortex","Motor_Cortex_Medial":"Motor_Cortex",
"Lateral_Motor_Cortex":"Motor_Cortex","Medial_Motor_Cortex":"Motor_Cortex","Primary_Motor_Cortex_L":"Motor_Cortex","Primary_Motor_Cortex_M":"Motor_Cortex",
"Frontal Cortex":"Frontal_Cortex","Cerebellum":"Cerebellum","Spinal_Cord_Cervical":"Cervical_Spinal_Cord","Cervical Spinal Cord":"Cervical_Spinal_Cord",
"Cervical_spinal_cord":"Cervical_Spinal_Cord","Spinal_cord_Cervical":"Cervical_Spinal_Cord","Lumbar Spinal Cord":"Lumbar_Spinal_Cord",
"Spinal_Cord_Lumbosacral":"Lumbar_Spinal_Cord","Lumbosacral_Spinal_Cord":"Lumbar_Spinal_Cord","Lumbar_spinal_cord":"Lumbar_Spinal_Cord"}
TISSUE_REMAP_NORM={norm(k):v for k,v in TISSUE_REMAP.items()}
for t in CANONICAL_TISSUES:TISSUE_REMAP_NORM[norm(t)]=t


def table(rows,path,columns=None):
    frame=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows,columns=columns)
    frame.to_csv(path,sep='\t',index=False,na_rep='NA')

def matrix(path,wanted):
    data={}
    with path.open() as f:
        subjects=f.readline().rstrip('\n').split('\t')[1:]
        if len(subjects)!=len(set(subjects)):raise ValueError(f'Duplicate subjects: {path}')
        for line in f:
            cells=line.rstrip('\n').split('\t');key=cells[0]
            if key not in wanted:continue
            if key in data:raise ValueError(f'Duplicate matrix row: {key}')
            if len(cells)-1!=len(subjects):raise ValueError(f'Wrong row length: {key}')
            data[key]=pd.to_numeric(pd.Series(cells[1:]),errors='coerce').to_numpy(float)
    return subjects,data

def ranked_pairs(rows):
    keep=(rows[['n_ref_ref','n_het','n_hom_alt']]>=5).all(axis=1)&np.isfinite(pd.to_numeric(rows.abs_delta_pp,errors="coerce"))
    rank=rows.loc[keep].sort_values(['abs_delta_pp','FDR','circ_id','snpid'],ascending=[False,True,True,True]).reset_index(drop=True)
    rank.insert(0,'rank',np.arange(1,len(rank)+1))
    return rank

def select_top(rank):
    return rank.loc[rank.groupby('circ_id',sort=False).cumcount()<2].head(10)

def sample_metadata(meta,tissue,subjects):
    patterns={'Motor_Cortex':'Motor Cortex Lateral|Motor Cortex Medial|Lateral Motor Cortex|Medial Motor Cortex|Primary Motor Cortex L|Primary Motor Cortex M|Cortex_Motor_Unspecified|Cortex_Motor_BA4|BA4 Motor Cortex|Lateral_motor_cortex|Motor Cortex|BA4',
              'Frontal_Cortex':'Frontal Cortex|Frontal','Cervical_Spinal_Cord':'Spinal_Cord_Cervical|Cervical Spinal Cord|Cervical_spinal_cord|Spinal_cord_Cervical|Cervical',
              'Lumbar_Spinal_Cord':'Lumbar Spinal Cord|Spinal_Cord_Lumbosacral|Lumbosacral_Spinal_Cord|Lumbar_spinal_cord|Lumbar|Lumbosacral','Cerebellum':'Cerebellum'}
    m=meta.copy()
    if tissue in ('Cortex','Spinal_Cord'):
        m=m[m.tissue.map(lambda v:tissue_equal(v,tissue))]
    else:
        pattern=patterns.get(tissue,re.escape(tissue))
        m=m[m.apply(lambda r:r.astype(str).str.contains(pattern,case=False,regex=True).any(),axis=1)]
    m=m[m.externalsubjectid.isin(subjects)].copy()
    m['RIN_score']=m.iloc[:,[16,17]].apply(pd.to_numeric,errors='coerce').fillna(0).max(axis=1)
    m=m[(pd.to_numeric(m.post_mortem_interval_in_hours,errors='coerce')<=40)&(m.RIN_score>=3)]
    return m.sort_values('RIN_score',ascending=False).drop_duplicates('externalsubjectid')

def read_metadata(path):
    return pd.read_csv(path,sep='\t' if path.suffix.lower() in ('.tsv','.txt') else ',',
                       dtype=str,keep_default_na=False,encoding='utf-8-sig')

def answer_sample_metadata(meta,tissue,subjects,subject_column=None,sample_column=None):
    m=meta.copy()
    m.columns=m.columns.str.strip()
    if m.empty:raise ValueError('Answer ALS metadata/sample map is empty')
    subject_column=subject_column or ('externalsubjectid' if 'externalsubjectid' in m else m.columns[0])
    sample_column=sample_column or ('externalsampleid' if 'externalsampleid' in m else subject_column)
    for col in (subject_column,sample_column):
        if col not in m:raise ValueError(f'Metadata column not found: {col}')
    print(f'  Answer ALS ID columns: subject={subject_column}, sample={sample_column}',flush=True)
    if 'tissue' in m:m=m[m.tissue.map(lambda x:tissue_equal(x,tissue))]
    selected=pd.DataFrame({'externalsubjectid':m[subject_column].astype(str).str.strip(),
                           'externalsampleid':m[sample_column].astype(str).str.strip()})
    selected=selected[selected.externalsubjectid.isin(subjects)].copy()
    if selected.empty:
        raise ValueError('No Answer ALS metadata IDs match matrix subjects. Supply --sample-map or explicit --subject-column/--sample-column.')
    if selected.externalsampleid.eq('').any():raise ValueError('Missing RNA sample ID in Answer ALS sample map')
    selected=selected.drop_duplicates()
    if selected.externalsubjectid.duplicated().any() or selected.externalsampleid.duplicated().any():
        raise ValueError('Answer ALS requires one RNA sample per matrix subject and one subject per RNA sample. Supply --sample-map with the samples used to build the circ matrix; do not guess replicate IDs.')
    # iPSC samples have no postmortem interval filter. Retain upstream QTL cohort.
    return selected

def quant_rows(path,wanted):
    # Exact circ ID only; no tolerance matching and no missing-to-zero conversion.
    result={};duplicate=set()
    if not path.exists():return result
    with path.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            try:
                key=f"chr{chr_key(row['chrom'])}:{int(float(row['start']))}-{int(float(row['end']))}:{row['strand'].strip()}"
                if key not in wanted:continue
                if key in result:duplicate.add(key);continue
                values={k:float(row[k]) for k in ['circ_percent','circ_reads','linear_reads']}
                if not all(np.isfinite(v) for v in values.values()):continue
                if not 0<=values['circ_percent']<=1 or min(values['circ_reads'],values['linear_reads'])<0:continue
                result[key]=values
            except (KeyError,ValueError):continue
    for key in duplicate:result.pop(key,None)
    return result

def render(pair,records,out,gtf,loc):
    event=pair['circ_id'];snp=pair['snpid'];coords=parse_event_coordinates(event)
    chrom=coords['chr'];start=coords['start'];end=coords['end'];strand=coords['strand']
    flank=min(20000,max(2000,math.ceil((end-start)*.2)));left=max(0,start-flank);right=end+flank
    nb=min(2500,max(50,math.ceil((right-left)/25)));width=(right-left)/nb
    x=np.linspace(left,right,nb,endpoint=False)+width/2
    # Never silently change the cohort used to rank the percentage difference.
    for row in records:
        p,m,ps,ms=bigwig_strands(row['sample_dir'],chrom,left,right,nb)
        if p is None or m is None or len(p)!=nb or len(m)!=nb or not np.isfinite(p).all() or not np.isfinite(m).all():
            return None,f"unusable_bigwig:{row['subject_id']}"
        row.update(plus_track=p,minus_track=m,plus_source=ps,minus_source=ms)
    groups={g:[r for r in records if r['genotype']==g] for g in range(3)}
    if any(len(g)<5 for g in groups.values()):raise ValueError('Plot cohort violates minimum N')
    exons=load_reference_exons(gtf,chrom,left,right)
    title=f"{tissue_display(pair['tissue'])} | cQTL | {event} | {snp} | REF={pair['ref']} ALT={pair['alt']}\nFDR={pair['FDR']:.4g} | ALT/ALT − REF/REF circ percentage = {pair['delta_alt_minus_ref_pp']:.3f} percentage points"
    try:
        fig,boundary_data=cq_boundary.make_figure(
            records,chrom,start,end,strand,exons,title,out/'.cqtl_boundary_cache',
            snp_pos=loc[1] if chr_key(loc[0])==chr_key(chrom) else None,
            reference_window=(left,right),event=event)
    except (OSError,ValueError,RuntimeError) as error:
        return None,'boundary_plot_failed:'+str(error)
    prefix=safe_name(f"{pair['tissue']}__{event}__{snp}")+'__'+hashlib.sha256(f'{event}|{snp}'.encode()).hexdigest()[:8]
    fig.savefig(out/(prefix+'.density.pdf'),bbox_inches='tight');fig.savefig(out/(prefix+'.density.png'),dpi=180,bbox_inches='tight');plt.close(fig)
    cq_boundary.export(boundary_data, out/prefix)
    export=[{k:str(v) if isinstance(v,Path) else v for k,v in r.items() if k not in ('plus_track','minus_track')} for r in records]
    table(export,out/(prefix+'.subjects.tsv'))
    return prefix,None

def aligned_records(event,subjects,samples,phenotype,g,ref,alt,counted):
    raw_missing=0;matrix_discordant=0
    records=[]
    for i,sub in enumerate(subjects):
        sm=samples.get(sub)
        if sm is None or not np.isfinite(g[i]) or not np.isfinite(phenotype[i]):continue
        quant=sm['quant'].get(event)
        if quant is None:raw_missing+=1;continue
        # Reject reconstructed sample mismatches rather than silently plotting a different RNA sample.
        if not np.isclose(quant['circ_percent'],phenotype[i],rtol=1e-5,atol=1e-6):matrix_discordant+=1;continue
        records.append(dict(subject_id=sub,sample_id=sm['sample_id'],sample_dir=sm['sample_dir'],primary_reads=sm['primary_reads'],genotype=int(g[i]),ref=ref,alt=alt,plink_counted_allele=counted,quant_file=str(sm['sample_dir']/'circularRNA_known_circ_percentage.txt'),**quant))
    return records,raw_missing,matrix_discordant

def main():
    ap=argparse.ArgumentParser(description='Top 10 significant cQTL circ percentage differences; max 2 SNPs per circ; N>=5 in all three genotype groups.')
    ap.add_argument('--experiment',choices=['target_ALS','answer_ALS'],default=os.environ.get('EXPERIMENT','target_ALS'))
    ap.add_argument('--root',type=Path)
    ap.add_argument('--metadata',type=Path)
    ap.add_argument('--sample-map',type=Path,help='Answer ALS TSV with externalsubjectid and externalsampleid; one sample per matrix subject.')
    ap.add_argument('--subject-column',help='Answer ALS metadata subject ID column.')
    ap.add_argument('--sample-column',help='Answer ALS metadata RNA sample ID column.')
    ap.add_argument('--processed-root',type=Path,help='Directory containing per-sample RNAseq folders; use with one tissue.')
    ap.add_argument('--outdir',type=Path)
    ap.add_argument('--tissues',nargs='+',default=None)
    ap.add_argument('--gtf',type=Path,default=Path.home()/'donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf')
    ap.add_argument('--no-plots',action='store_true')
    a=ap.parse_args()
    global EXPERIMENT,PROCESSED_ROOT
    EXPERIMENT=a.experiment
    if EXPERIMENT not in ('target_ALS','answer_ALS'):ap.error('Invalid EXPERIMENT environment variable')
    a.root=(a.root or Path.home()/'donglab/data'/EXPERIMENT).expanduser()
    defaults=['iPSC_Motor_Neuron'] if EXPERIMENT=='answer_ALS' else CANONICAL_TISSUES[:5]
    a.tissues=a.tissues or os.environ.get('QTL_TISSUES','').split() or defaults
    if a.processed_root and len(a.tissues)!=1:ap.error('--processed-root requires exactly one tissue')
    PROCESSED_ROOT=a.processed_root.expanduser() if a.processed_root else None
    if PROCESSED_ROOT and not PROCESSED_ROOT.is_dir():ap.error(f'Processed root not found: {PROCESSED_ROOT}')
    if EXPERIMENT=='target_ALS' and any([a.sample_map,a.subject_column,a.sample_column]):
        ap.error('--sample-map/--subject-column/--sample-column are for --experiment answer_ALS')
    if a.sample_map and any([a.metadata,a.subject_column,a.sample_column]):
        ap.error('Use --sample-map by itself, or --metadata with optional ID columns')
    a.outdir=(a.outdir or a.root/'QTL/cqtl_circ_difference').expanduser()
    a.outdir.mkdir(parents=True,exist_ok=True)
    print(f'Experiment: {EXPERIMENT} | root: {a.root} | tissues: {a.tissues}',flush=True)
    print('Loading gene symbols from',a.gtf,flush=True)
    gene_index=load_gene_index(a.gtf);gene_annotation_cache={}
    metadata_path=(a.sample_map or a.metadata or a.root/('aals_dataportal_datatable.csv' if EXPERIMENT=='answer_ALS' else 'targetALS_rnaseq_metadata.csv')).expanduser()
    print('Metadata:',metadata_path,flush=True)
    meta=read_metadata(metadata_path) if EXPERIMENT=='answer_ALS' else pd.read_csv(metadata_path,low_memory=False)
    if a.sample_map and not {'externalsubjectid','externalsampleid'}.issubset(meta.columns):
        ap.error('--sample-map requires externalsubjectid and externalsampleid columns')
    allele_map=defaultdict(set)
    with (a.root/'QTL/plink/joint_all_chrs_matrixEQTL.raw').open() as f:
        for token in f.readline().split()[6:]:
            m=re.fullmatch(r'(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)',token,re.I)
            if m:
                c,p,ref,alt,counted=m.groups();allele_map[(chr_key(c),int(p))].add((ref.upper(),alt.upper(),counted.upper()))
    tops=[];audit=[];jobs=[]
    for tissue in a.tissues:
        print('START',tissue,flush=True)
        q=a.root/tissue/'cQTL';out=a.outdir/tissue;out.mkdir(parents=True,exist_ok=True)
        source=q/'results'/f'{tissue}_cQTL.FDR0.05.txt'
        sig=pd.read_csv(source,sep='\t',dtype={'circ_id':str,'snpid':str})
        sig['FDR']=pd.to_numeric(sig.FDR,errors='coerce')
        sig=sig[sig.FDR.between(0,.05)].sort_values('FDR').drop_duplicates(['circ_id','snpid'])
        subjects,phenos=matrix(q/f'circ_{tissue}.txt',set(sig.circ_id))
        gs,genotypes=matrix(q/f'snp_{tissue}.txt',set(sig.snpid))
        if set(gs)!=set(subjects):raise ValueError(f'{tissue}: phenotype/genotype subjects differ')
        order=[gs.index(s) for s in subjects]
        loc={}
        with (q/'snp_location.txt').open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                key=r['snpid'];v=(chr_key(r['chr']),int(r['pos']))
                if key in loc and loc[key]!=v:raise ValueError(f'Ambiguous location: {key}')
                loc[key]=v
        mt=answer_sample_metadata(meta,tissue,subjects,a.subject_column,a.sample_column) if EXPERIMENT=='answer_ALS' else sample_metadata(meta,tissue,subjects)
        samples={};sample_audit=[];wanted=set(sig.circ_id)
        for row in mt.to_dict('records'):
            sid=str(row['externalsampleid']);sub=str(row['externalsubjectid']);sd=resolve_sample_dir(a.root,tissue,sid)
            reason=''
            if sd is None:reason='missing_sample_directory'
            elif not all((sd/name).is_file() for name in ['STAR.Aligned.sortedByCoord.out.plus.normalized.bw','STAR.Aligned.sortedByCoord.out.minus.normalized.bw']):reason='missing_bigwig'
            n=primary_read_count(sd) if not reason else None
            if not reason and not n:reason='missing_primary_read_count'
            if reason:sample_audit.append(dict(subject_id=sub,sample_id=sid,reason=reason));continue
            quant=quant_rows(sd/'circularRNA_known_circ_percentage.txt',wanted)
            samples[sub]=dict(sample_id=sid,sample_dir=sd,primary_reads=n,quant=quant)
        for sub in set(subjects)-set(mt.externalsubjectid):sample_audit.append(dict(subject_id=sub,sample_id='',reason='missing_eligible_metadata'))
        table(sample_audit,out/'excluded_samples.tsv',columns=['subject_id','sample_id','reason'])
        print(f'  Matrix subjects: {len(subjects)} | selected metadata: {len(mt)} | usable samples: {len(samples)}',flush=True)
        if not samples:raise ValueError(f'{tissue}: no usable RNA samples. See {out}/excluded_samples.tsv; check --processed-root, ID mapping, BAM and BigWig filenames.')
        rows=[];excluded=[];raw_missing=0;matrix_discordant=0
        for pair in sig.to_dict('records'):
            event=pair['circ_id'];snp=pair['snpid'];aa=allele_map.get(loc.get(snp),set())
            reason=None
            if event not in phenos or snp not in genotypes:reason='missing_matrix_row'
            elif parse_event_coordinates(event) is None:reason='unparseable_circ_coordinates'
            elif len(aa)!=1:reason='missing_or_ambiguous_alleles'
            if reason:excluded.append(dict(pair,reason=reason));continue
            ref,alt,counted=next(iter(aa))
            if ref==alt or counted not in (ref,alt):excluded.append(dict(pair,reason='unsupported_alleles'));continue
            g=genotypes[snp][order]
            if np.any(np.isfinite(g)&~np.isin(g,[0,1,2])):raise ValueError(f'Non-hardcall genotype: {snp}')
            if counted==ref:g=2-g
            records,nmissing,ndiscordant=aligned_records(event,subjects,samples,phenos[event],g,ref,alt,counted)
            raw_missing+=nmissing;matrix_discordant+=ndiscordant
            result=dict(pair,tissue=tissue,ref=ref,alt=alt,plink_counted_allele=counted)
            result.update(annotate_circ_gene(event,gene_index,gene_annotation_cache))
            for value,label in [(0,'ref_ref'),(1,'het'),(2,'hom_alt')]:
                v=[r['circ_percent'] for r in records if r['genotype']==value]
                result['n_'+label]=len(v);result['mean_'+label+'_pct']=100*np.mean(v) if v else np.nan
            result['delta_alt_minus_ref_pp']=result['mean_hom_alt_pct']-result['mean_ref_ref_pct'];result['abs_delta_pp']=abs(result['delta_alt_minus_ref_pp'])
            rows.append(result)
        columns=list(dict.fromkeys(list(sig.columns)+['gene_name','gene_id','tissue','ref','alt','plink_counted_allele','n_ref_ref','n_het','n_hom_alt','mean_ref_ref_pct','mean_het_pct','mean_hom_alt_pct','delta_alt_minus_ref_pp','abs_delta_pp']))
        allrows=pd.DataFrame(rows,columns=columns);rank=ranked_pairs(allrows)
        table(rank,out/'ranked_qualifying.tsv');table(rank.drop_duplicates('circ_id'),out/'best_snp_per_circ.tsv')
        rejected=allrows[(allrows[['n_ref_ref','n_het','n_hom_alt']]<5).any(axis=1)].copy();rejected['reason']='fewer_than_5_usable_subjects_in_at_least_one_genotype'
        table(rejected,out/'excluded_genotype_counts.tsv');table(excluded,out/'excluded_mapping.tsv',columns=list(sig.columns)+['reason'])
        selected=[];counts=defaultdict(int);plot_exclusions=[]
        for pair in rank.to_dict('records'):
            if counts[pair['circ_id']]>=2:continue
            prefix=None
            if not a.no_plots:
                g=genotypes[pair['snpid']][order]
                if pair['plink_counted_allele']==pair['ref']:g=2-g
                records,_,_=aligned_records(pair['circ_id'],subjects,samples,phenos[pair['circ_id']],g,pair['ref'],pair['alt'],pair['plink_counted_allele'])
                prefix,reason=render(pair,records,out,a.gtf,loc[pair['snpid']])
                if reason:plot_exclusions.append(dict(pair,reason=reason));continue
            selected.append(pair);counts[pair['circ_id']]+=1
            jobs.append(dict(pair,prefix=prefix,outdir=str(out),status='ranking_only' if a.no_plots else 'plotted'))
            if len(selected)==10:break
        table(plot_exclusions,out/'excluded_plot_inputs.tsv',columns=['circ_id','snpid','reason'])
        pd.DataFrame(selected,columns=columns)[['circ_id','snpid']].to_csv(out/'top_pairs.tsv',sep='\t',header=False,index=False)
        tops.extend(selected);audit.append(dict(tissue=tissue,significant_pairs=len(sig),qualifying_pairs=len(rank),selected_for_plots=len(selected),plots_written=0 if a.no_plots else len(selected),mapping_excluded=len(excluded),excluded_genotype_counts=len(rejected),missing_raw_quant_pair_subjects=raw_missing,discordant_matrix_pair_subjects=matrix_discordant,plot_input_exclusions=len(plot_exclusions)))
        print(tissue,'selected',len(selected),'qualifying',len(rank),flush=True)
    top_columns=list(dict.fromkeys(['rank','tissue','circ_id','gene_name','gene_id','snpid']+columns))
    table(tops,a.outdir/'top10_per_tissue.tsv',columns=top_columns)
    table(audit,a.outdir/'audit.tsv');(a.outdir/'plot_jobs.json').write_text(json.dumps(jobs,indent=2,default=str))
    print('COMPLETE',a.outdir,flush=True)

if __name__=='__main__':main()
PYTHON
