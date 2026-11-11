#!/bin/bash
# BOXPLOT_VERSION: EXACT_COPY_OF_ORIGINAL_BASE_R_PARAMETERS__MATCHED_TRACK_SUBJECTS_ONLY
# GENCODE_TRACK_VERSION: ALL_TRANSCRIPTS_TRANSCRIPT_RESOLVED
#SBATCH --job-name=run_coloc_plot
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/run_coloc_plot.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G

TYPE="${1:-}"
[[ "$TYPE" == "circ eQTL" ]] && TYPE="circ_eQTL"
ROOT="/home/zw529/donglab/data/target_ALS"; MR_DIR="${ROOT}/MR"
METADATA="${ROOT}/targetALS_rnaseq_metadata.csv"
GTF="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Annotation/gencode/gencode.v49.annotation.gtf"
RAW_FILE="${ROOT}/QTL/plink/joint_all_chrs_matrixEQTL.raw"
FLANK_FRAC=0.20; FLANK_MIN=2000; FLANK_MAX=20000; BIGWIG_BIN_BP=25; MAX_TRACK_POINTS=2500
CIRC_COORD_TOL=3; MIN_MEAN_JUNCTION_RPM=0.01; MAX_JUNCTIONS=30; MAX_HITS=0; TISSUE_FILTER=""; OVERWRITE=0; DPI=180

set -euo pipefail
case "$TYPE" in eQTL|sQTL|cQTL|circ_eQTL) ;; *) echo "Usage: sbatch $0 {eQTL|sQTL|cQTL|circ_eQTL}"; exit 1 ;; esac
COLOC="${MR_DIR}/${TYPE}_SuSiE_coloc_summary.tsv"; OUTDIR="${MR_DIR}/coloc_raw_plots/${TYPE}"; mkdir -p "$OUTDIR"
for f in "$COLOC" "$METADATA" "$GTF" "$RAW_FILE"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done
module load deepTools 2>/dev/null || { echo "ERROR: could not load deepTools"; exit 1; }
PYTHON="$(command -v python)"; command -v samtools >/dev/null || { echo "ERROR: samtools unavailable"; exit 1; }
"$PYTHON" -c 'import numpy,matplotlib,pyBigWig,scipy' || { echo "ERROR: Python plotting stack unavailable"; exit 1; }
export TYPE ROOT METADATA COLOC OUTDIR GTF RAW_FILE FLANK_FRAC FLANK_MIN FLANK_MAX BIGWIG_BIN_BP MAX_TRACK_POINTS CIRC_COORD_TOL
export MIN_MEAN_JUNCTION_RPM MAX_JUNCTIONS MAX_HITS TISSUE_FILTER OVERWRITE DPI
"$PYTHON" - <<'PY'

import os,re,csv,sys,math,warnings,subprocess,hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle,PathPatch
from matplotlib.path import Path as MplPath
import pyBigWig
warnings.filterwarnings("ignore")

GENOTYPE_LABELS={0:"Ref/Ref",1:"Het",2:"Hom Alt"}
TRACK_COLORS={0:"#228B22",1:"#4682B4",2:"#C75B7A"}
CANONICAL_TISSUES=["Cerebellum","Frontal_Cortex","Cervical_Spinal_Cord","Lumbar_Spinal_Cord","Motor_Cortex"]
if os.environ.get("QTL_TISSUES"):CANONICAL_TISSUES=os.environ["QTL_TISSUES"].split()

def norm(s):return re.sub(r"[^a-z0-9]","",str(s).lower())
def chr_key(s):return re.sub(r"^chr","",str(s),flags=re.I).upper()
def safe_name(s):return re.sub(r"[^A-Za-z0-9._+-]+","_",str(s))[:180]
def strip_gene_version(x):return re.sub(r"\.\d+$","",str(x))
def tissue_display(x):return str(x).replace("_"," ")
def missing(v):return v is None or str(v).strip().lower() in {"","na","nan","none",".","null"}
def truthy(v):return str(v).strip().lower() in {"true","t","1","yes","y","pass","passed"}
def detect_delimiter(line):return "\t" if "\t" in line else None
def split_line(line,d):return line.rstrip("\r\n").split("\t") if d=="\t" else line.split()
def as_float(x):
    try:v=float(x);return v if np.isfinite(v) else None
    except:return None
def format_p(x):return "NA" if as_float(x) is None else f"{float(x):.3e}"
def format_max_sig(v,sig=5):
    if v==0:return "0"
    d=max(0,sig-int(math.floor(math.log10(abs(v))))-1);return f"{v:.{d}f}".rstrip("0").rstrip(".")
def genomic_formatter(x,pos=None):return f"{format_max_sig(x/1e6)} Mb" if abs(x)>=1e6 else f"{format_max_sig(x/1e3)} kb" if abs(x)>=1e3 else format_max_sig(x)
def density_formatter(y,pos=None):return "0" if abs(y)<1e-12 else format_max_sig(abs(y),4)
def find_header(headers,aliases,required=True):
    lookup={norm(h):h for h in headers}
    for a in aliases:
        if norm(a) in lookup:return lookup[norm(a)]
    for h in headers:
        if any(norm(a) in norm(h) for a in aliases if norm(a)):return h
    if required:raise KeyError(f"Could not find {aliases}. Available: {headers}")
    return None

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
def canonical_tissue(x):x=str(x).strip();return TISSUE_REMAP_NORM.get(norm(x),x)
def tissue_equal(a,b):
    if b=="Cortex":return bool(re.search(r"(?:^|[\W_])cortex(?:$|[\W_])",str(a),re.I))
    if b=="Spinal_Cord":return bool(re.search(r"(?:^|[\W_])spinal[\W_]*cord(?:$|[\W_])",str(a),re.I))
    return canonical_tissue(a)==canonical_tissue(b)

def extract_matrix_rows(path,wanted):
    wanted=set(map(str,wanted));wanted_nv={strip_gene_version(x) for x in wanted};out={}
    with open(path) as f:
        first=f.readline();d=detect_delimiter(first);samples=split_line(first,d)[1:]
        for line in f:
            if not line.strip():continue
            x=split_line(line,d)
            if len(x)<2:continue
            rid=x[0].strip()
            if rid not in wanted and strip_gene_version(rid) not in wanted_nv:continue
            vals=x[1:]+[""]*max(0,len(samples)-len(x[1:]));out[rid]=dict(zip(samples,vals[:len(samples)]))
    return out
def find_matrix_row(rows,event):
    if event in rows:return rows[event]
    e=strip_gene_version(event);return next((v for k,v in rows.items() if strip_gene_version(k)==e),None)

def load_covariate_status(path):
    if not path.exists():return {}
    with open(path) as f:
        first=f.readline();d=detect_delimiter(first);samples=split_line(first,d)[1:]
        for line in f:
            x=split_line(line,d)
            if x and x[0].strip()=="is_als":
                vals=x[1:]+[""]*max(0,len(samples)-len(x[1:]))
                out={}
                for s,v in zip(samples,vals[:len(samples)]):
                    try:out[s]=1 if float(v)>.5 else 0
                    except:out[s]=0
                return out
    return {}

def load_location_table(path,kind,TYPE):
    with open(path) as f:
        first=f.readline();d=detect_delimiter(first);fields=split_line(first,d)
        is_data=len(fields)>=3 and re.match(r"^(chr)?[0-9XYM]+$",fields[1],re.I) and re.match(r"^\d+$",fields[2])
        if is_data:
            headers=["id","chr","start"]+(["end"] if len(fields)>=4 else [])+(["strand"] if len(fields)>=5 else [])
            while len(headers)<len(fields):headers.append(f"extra{len(headers)}")
            lines=[fields]
        else:headers,lines=fields,[]
        lines += [split_line(x,d) for x in f if x.strip()]
    rows=[]
    for x in lines:x+=[""]*max(0,len(headers)-len(x));rows.append(dict(zip(headers,x)))
    if kind=="snp":
        i=find_header(headers,["snpid","snp_id","rsid","id"]);c=find_header(headers,["chr","chrom","chromosome"]);p=find_header(headers,["pos","position","start"])
        return {str(r[i]):(str(r[c]),int(float(r[p]))) for r in rows if not any(missing(r.get(z)) for z in [i,c,p])}
    aliases={"eQTL":["geneid","gene_id","gene","id"],"sQTL":["junction_id","junctionid","splicing_id","event_id","phenotype_id","id"],"cQTL":["circ_id","circid","circrna_id","event_id","phenotype_id","id"]}
    i=find_header(headers,aliases["eQTL" if TYPE=="circ_eQTL" else TYPE]);c=find_header(headers,["chr","chrom","chromosome"]);s=find_header(headers,["start","left","begin","donor"])
    e=find_header(headers,["end","right","stop","acceptor"]);st=find_header(headers,["strand"],False);out={}
    for r in rows:
        try:out[str(r[i])]={"chr":str(r[c]),"start":int(float(r[s])),"end":int(float(r[e])),"strand":str(r[st]) if st and not missing(r.get(st)) else None}
        except:pass
    return out
def parse_event_coordinates(event):
    for p in [r"(chr[^:]+):(\d+)-(\d+):([+-])",r"(chr[^:]+):(\d+):(\d+):([+-])",r"(chr[^:]+):(\d+)-(\d+)",r"(chr[^:]+):(\d+):(\d+)"]:
        m=re.search(p,str(event))
        if m:return {"chr":m.group(1),"start":int(m.group(2)),"end":int(m.group(3)),"strand":m.group(4) if len(m.groups())>=4 else None}
    return None
def find_location(loc,event):
    if event in loc:return loc[event]
    e=strip_gene_version(event);return next((v for k,v in loc.items() if strip_gene_version(k)==e),parse_event_coordinates(event))

def make_matched_boxplot(complete,out,event,snp,symbol,allele,TYPE):
    if len(complete)<5:return False
    import tempfile,shlex
    with tempfile.TemporaryDirectory(prefix="matched_boxplot_") as td:
        td=Path(td); data=td/"data.tsv"; rfile=td/"boxplot.R"
        with open(data,"w",newline="") as fh:
            w=csv.writer(fh,delimiter="\t");w.writerow(["phenotype","SNP","is_als"])
            for r in complete:w.writerow([r["phenotype"],r["genotype"],r.get("is_als",0)])
        rcode=r'''suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(TRUE)
d <- fread(args[1], header=TRUE)
out <- args[2]; event <- args[3]; snp <- args[4]; symbol <- args[5]; ref <- args[6]; alt <- args[7]; qtype <- args[8]
df <- data.frame(
  phenotype=as.numeric(d$phenotype),
  SNP=as.numeric(d$SNP),
  p_col=ifelse(as.numeric(d$is_als) > 0.5, "red", "#9932CC"),
  p_pch=ifelse(as.numeric(d$is_als) > 0.5, 16, 18),
  p_cex=ifelse(as.numeric(d$is_als) > 0.5, 0.7, 1.1),
  stringsAsFactors=FALSE
)
df <- df[!is.na(df$SNP) & !is.na(df$phenotype), ]

get_p <- function(formula, data) {
  tryCatch({
    fit <- lm(formula, data=data)
    formatC(summary(fit)$coefficients[2,4], format="e", digits=3)
  }, error=function(e) "NA")
}

ylab <- if (qtype=="circ_eQTL") "circRNA abundance (rank inverse-normal RPM)" else if (qtype=="eQTL") "Expression (Z-score)" else if (qtype=="sQTL") "Splicing Level (PSI)" else "circRNA Percentage (%)"

pdf(out, width=12, height=5)
par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))

add_counts <- table(factor(df$SNP, levels=0:2))
add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)

boxplot(phenotype ~ SNP_f, data=df, col="lightgreen", outline=FALSE,
        ylab=ylab, xlab="Genotype",
        main=paste0("Additive Model\n(p = ", get_p(phenotype ~ SNP, df), ")"))

points(jitter(as.numeric(df$SNP_f), amount=0.15), df$phenotype,
       pch=df$p_pch, col=df$p_col, cex=df$p_cex)

legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
       pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

df$dom_val <- ifelse(df$SNP > 0, 1, 0)
dom_counts <- table(factor(df$dom_val, levels=0:1))
dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)

boxplot(phenotype ~ dom_f, data=df, col="lightblue", outline=FALSE,
        ylab=ylab, xlab="Genotype",
        main=paste0("Dominant Model\n(p = ", get_p(phenotype ~ dom_val, df), ")"))

points(jitter(as.numeric(df$dom_f), amount=0.15), df$phenotype,
       pch=df$p_pch, col=df$p_col, cex=df$p_cex)

legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
       pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

df$rec_val <- ifelse(df$SNP == 2, 1, 0)
rec_counts <- table(factor(df$rec_val, levels=0:1))
rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)

boxplot(phenotype ~ rec_f, data=df, col="lightpink", outline=FALSE,
        ylab=ylab, xlab="Genotype",
        main=paste0("Recessive Model\n(p = ", get_p(phenotype ~ rec_val, df), ")"))

points(jitter(as.numeric(df$rec_f), amount=0.15), df$phenotype,
       pch=df$p_pch, col=df$p_col, cex=df$p_cex)

legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
       pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

if (qtype=="eQTL") {
  mtext(paste("Gene:", symbol, "(", event, ") | SNP:", snp, "| REF:", ref, "| ALT:", alt),
        outer=TRUE, cex=1.0, font=2, line=0.5)
} else if (qtype=="sQTL") {
  mtext(paste("Junction:", event, "| SNP:", snp, "| REF:", ref, "| ALT:", alt),
        outer=TRUE, cex=1.1, font=2, line=0.5)
} else {
  mtext(paste("circRNA:", event, "| SNP:", snp, "| REF:", ref, "| ALT:", alt),
        outer=TRUE, cex=1.0, font=2, line=0.5)
}
dev.off()
'''
        rfile.write_text(rcode)
        args=[str(data),str(out),str(event),str(snp),str(symbol),str(allele["ref"]),str(allele["alt"]),TYPE]
        cmd="module load R >/dev/null 2>&1; Rscript "+shlex.quote(str(rfile))+" "+" ".join(shlex.quote(x) for x in args)
        try:
            subprocess.run(["bash","-lc",cmd],check=True)
        except subprocess.CalledProcessError as e:
            print(f"  BOXPLOT FAILED: Rscript exited {e.returncode}");return False
    print(f"  BOXPLOT REGENERATED (exact original R aesthetics): N={len(complete)} -> {out}")
    return True

processed_root_cache={};processed_sample_cache={};primary_read_cache={}
def processed_roots_for_tissue(ROOT,tissue):
    t=canonical_tissue(tissue)
    if t in processed_root_cache:return processed_root_cache[t]
    roots=[];preferred=ROOT/t/"RNAseq"/"Processed"
    if preferred.is_dir():roots.append(preferred)
    for p in ROOT.iterdir():
        if p.is_dir() and tissue_equal(p.name,t):
            q=p/"RNAseq"/"Processed"
            if q.is_dir() and q not in roots:roots.append(q)
    processed_root_cache[t]=roots;print(f"  RNAseq roots for {t}: "+(" | ".join(map(str,roots)) if roots else "NONE"));return roots
def resolve_sample_dir(ROOT,tissue,sid):
    key=(canonical_tissue(tissue),norm(sid))
    if key in processed_sample_cache:return processed_sample_cache[key]
    for root in processed_roots_for_tissue(ROOT,tissue):
        for x in [sid,sid.replace("-","_"),sid.replace("_","-")]:
            p=root/x
            if p.is_dir():processed_sample_cache[key]=p;return p
    for root in processed_roots_for_tissue(ROOT,tissue):
        rk=str(root)
        if rk not in processed_sample_cache:processed_sample_cache[rk]={norm(p.name):p for p in root.iterdir() if p.is_dir()}
        p=processed_sample_cache[rk].get(norm(sid))
        if p is not None:processed_sample_cache[key]=p;return p
    processed_sample_cache[key]=None;return None
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

leafcutter_cache={}
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

def genotype_panel_title(group,genotype,qtype):
    title=f"{GENOTYPE_LABELS[genotype]} (N={len(group)})"
    if qtype=="circ_eQTL":
        values=[r["circ_rpm"] for r in group if r["circ_rpm"] is not None]
        return title+(f" | mean BSJ RPM={np.mean(values):.4g}" if values else " | mean BSJ RPM=NA")
    if qtype!="cQTL":return title
    # Match run_cQTL_circ_difference.sh: average per-subject fractions,
    # then multiply by 100 for display; this is not a ratio of pooled reads.
    values=[r["circ_percent_file"] for r in group
            if r["circ_percent_file"] is not None
            and math.isfinite(r["circ_percent_file"]) and 0<=r["circ_percent_file"]<=1]
    title+=f" | mean circ={100*sum(values)/len(values):.3f}%" if values else " | mean circ=NA"
    if len(values)!=len(group):title+=f" (quant N={len(values)}/{len(group)})"
    return title

def circ_quant_file(sd):return sd/"circularRNA_known_circ_percentage.txt"
def read_circ_quant(path,chrom,start,end,strand,tol):
    if not path.exists():return None
    exact=None;nearby=[]
    with open(path) as f:
        for r in csv.DictReader(f,delimiter="\t"):
            try:
                s,e=int(float(r["start"])),int(float(r["end"]));st=r["strand"].strip()
                if chr_key(r["chrom"])!=chr_key(chrom) or (strand in {"+","-"} and st!=strand):continue
                q={"circ_reads":float(r["circ_reads"]),"linear_reads":float(r["linear_reads"]),"circ_percent":float(r["circ_percent"]),"matched_start":s,"matched_end":e,"match_type":None}
                if s==start and e==end:q["match_type"]="exact";exact=q;break
                if abs(s-start)<=tol and abs(e-end)<=tol:q["match_type"]="tolerance";nearby.append(q)
            except:continue
    if exact:return exact
    if not nearby:return None
    nearby.sort(key=lambda r:abs(r["matched_start"]-start)+abs(r["matched_end"]-end));return nearby[0]

def parse_gtf_attrs(s):
    out={}
    for x in s.strip().split(";"):
        p=x.strip().split(" ",1)
        if len(p)==2:out[p[0]]=p[1].strip().strip('"')
    return out
def load_reference_exons(gtf,chrom,start,end):
    out=[]
    with open(gtf) as f:
        for l in f:
            if l.startswith("#"):continue
            x=l.rstrip().split("\t")
            if len(x)<9 or x[2]!="exon" or chr_key(x[0])!=chr_key(chrom):continue
            s,e=int(x[3])-1,int(x[4])
            if e<start or s>end:continue
            a=parse_gtf_attrs(x[8])
            tid=a.get("transcript_id","")
            if not tid:continue
            out.append({
                "start":s,
                "end":e,
                "strand":x[6],
                "gene":a.get("gene_name",a.get("gene_id","Unknown")),
                "gene_id":a.get("gene_id",""),
                "transcript_id":tid,
                "transcript_name":a.get("transcript_name","")
            })
    return out

def draw_reference_track(ax,exons,start,end):
    transcripts=defaultdict(list)
    info={}
    for e in exons:
        tid=e["transcript_id"]
        transcripts[tid].append((e["start"],e["end"]))
        info[tid]={
            "gene":e["gene"],
            "strand":e["strand"],
            "transcript_name":e["transcript_name"]
        }

    if not transcripts:
        ax.text(.5,.45,"No GENCODE v49 transcripts",transform=ax.transAxes,ha="center",fontsize=9)
        ax.set_ylim(-.45,.55)
    else:
        order=sorted(
            transcripts,
            key=lambda tid:(
                min(s for s,e in transcripts[tid]),
                info[tid]["gene"],
                tid
            )
        )

        for y,tid in enumerate(order):
            ints=sorted(transcripts[tid])
            tx_start=max(start,min(s for s,e in ints))
            tx_end=min(end,max(e for s,e in ints))

            # Transcript backbone.
            ax.plot([tx_start,tx_end],[y,y],color="black",lw=.7)

            # Preserve every exon from this specific transcript.
            for s,e in ints:
                s=max(s,start);e=min(e,end)
                if e>s:
                    ax.add_patch(Rectangle((s,y-.14),e-s,.28,facecolor="black",edgecolor="black",lw=.35))

            gene=info[tid]["gene"]
            strand=info[tid]["strand"]
            tname=info[tid]["transcript_name"]

            if tname and tname!=tid:
                label=f"{gene} | {tname} | {tid} ({strand})"
            else:
                label=f"{gene} | {tid} ({strand})"

            ax.text(start,y+.19,label,fontsize=6.3,ha="left",va="bottom")

        ax.set_ylim(-.45,len(order)-.05)

    ax.set_xlim(start,end)
    ax.set_yticks([])
    ax.set_ylabel("GENCODE v49",fontsize=12,labelpad=13)
    ax.xaxis.set_major_formatter(FuncFormatter(genomic_formatter))
    ax.locator_params(axis="x",nbins=5)
    ax.tick_params(axis="x",labelsize=11)
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

def genotype_group(v,flip=False):
    try:x=float(v)
    except:return None
    if flip:x=2-x
    return int(np.clip(np.rint(x),0,2)) if np.isfinite(x) and -.1<=x<=2.1 else None

TYPE,ROOT=os.environ["TYPE"],Path(os.environ["ROOT"])
METADATA,COLOC,OUTDIR,GTF,RAW_FILE=map(Path,[os.environ["METADATA"],os.environ["COLOC"],os.environ["OUTDIR"],os.environ["GTF"],os.environ["RAW_FILE"]])
FLANK_FRAC=float(os.environ["FLANK_FRAC"]);FLANK_MIN=int(os.environ["FLANK_MIN"]);FLANK_MAX=int(os.environ["FLANK_MAX"]);BIGWIG_BIN_BP=int(os.environ["BIGWIG_BIN_BP"]);MAX_TRACK_POINTS=int(os.environ["MAX_TRACK_POINTS"])
CIRC_COORD_TOL=int(os.environ["CIRC_COORD_TOL"]);MIN_MEAN_JUNCTION_RPM=float(os.environ["MIN_MEAN_JUNCTION_RPM"]);MAX_JUNCTIONS=int(os.environ["MAX_JUNCTIONS"]);MAX_HITS=int(os.environ["MAX_HITS"])
TISSUE_FILTER=os.environ["TISSUE_FILTER"].strip();OVERWRITE=int(os.environ["OVERWRITE"]);DPI=int(os.environ["DPI"])

raw_alleles={}
with open(RAW_FILE) as f:
    for col in f.readline().split():
        m=re.match(r"^(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)_([ACGT]+)$",col,re.I)
        if m:
            c,p,r,a,counted=m.groups();raw_alleles[(chr_key(c),int(p))]={"ref":r.upper(),"alt":a.upper(),"counted":counted.upper()}
print(f"PLINK RAW allele mappings: {len(raw_alleles):,}")

with open(COLOC,newline="") as f:
    r=csv.DictReader(f,delimiter="\t");summary_headers=r.fieldnames or []
    if TYPE=="circ_eQTL":
        if "max_PP_H4" not in summary_headers:sys.exit("ERROR: circ-eQTL summary lacks max_PP_H4")
        hits=[x for x in r if as_float(x.get("max_PP_H4")) is not None and float(x["max_PP_H4"])>0.8]
    else:
        if "H4_reference_pass" not in summary_headers:sys.exit("ERROR: summary lacks H4_reference_pass")
        hits=[x for x in r if truthy(x.get("H4_reference_pass"))]
tissue_col=find_header(summary_headers,["tissue"])
hits=[h for h in hits if canonical_tissue(h.get(tissue_col,"")) in CANONICAL_TISSUES]
if TISSUE_FILTER:hits=[h for h in hits if tissue_equal(h.get(tissue_col,""),TISSUE_FILTER)]
if MAX_HITS>0:hits=hits[:MAX_HITS]
if not hits:sys.exit("No passing coloc hits to plot.")
if TYPE=="eQTL":event_col=find_header(summary_headers,["geneid","gene_id","gene"]);symbol_col=find_header(summary_headers,["gene_symbol","symbol"],False)
elif TYPE=="sQTL":event_col=find_header(summary_headers,["junction_id","junctionid","splicing_id","splice_id","event_id","phenotype_id"]);symbol_col=find_header(summary_headers,["gene_symbol","symbol","gene"],False)
elif TYPE=="circ_eQTL":event_col=find_header(summary_headers,["circ_id"]);symbol_col=find_header(summary_headers,["host_gene","gene_symbol","symbol"],False)
else:event_col=find_header(summary_headers,["circ_id","circid","circrna_id","event_id","phenotype_id"]);symbol_col=find_header(summary_headers,["gene_symbol","symbol","gene"],False)
print("="*60);print(f"QTL TYPE        : {TYPE}");print(f"H4 PASSING HITS : {len(hits)}");print("FILTER          : "+("max_PP_H4 > 0.8" if TYPE=="circ_eQTL" else "H4_reference_pass == TRUE"));print("="*60)

with open(METADATA,newline="") as f:r=csv.DictReader(f);meta_headers,metadata=r.fieldnames or [],list(r)
subject_col=find_header(meta_headers,["externalsubjectid","subject_id","subjectid"]);sample_col=find_header(meta_headers,["externalsampleid","sample_id","sampleid"])
meta_tissue_col=find_header(meta_headers,["tissue"],False);rin_col=find_header(meta_headers,["rin","rna_integrity_number"],False)
def get_rin(r):
    try:return float(r.get(rin_col,"")) if rin_col else -1
    except:return -1
def metadata_for_qtl_id(qid,tissue):
    rows=metadata;target=canonical_tissue(tissue)
    if meta_tissue_col:
        z=[r for r in rows if tissue_equal(r.get(meta_tissue_col,""),target)]
        if z:rows=z
    z=[r for r in rows if str(r.get(sample_col,"")).strip()==str(qid)]
    if z:return max(z,key=get_rin)
    z=[r for r in rows if str(r.get(subject_col,"")).strip()==str(qid)]
    return max(z,key=get_rin) if z else None

def possible_snps(hit):
    c=[]
    for col in ["best_candidate_SNP","best_candidate_snp"]:
        if col in hit and not missing(hit[col]):c.append(hit[col])
    for col in {"circ_eQTL":["top_eQTL_SNP","top_QTL_SNP"],"eQTL":["top_eQTL_SNP","top_QTL_SNP"],"sQTL":["top_sQTL_SNP","top_QTL_SNP"],"cQTL":["top_cQTL_SNP","top_QTL_SNP"]}[TYPE]:
        if col in hit and not missing(hit[col]):c.append(hit[col])
    for k,v in hit.items():
        if "topsnp" in norm(k) and "gwas" not in norm(k) and not missing(v):c.append(v)
    if "candidate_SNPs" in hit and not missing(hit["candidate_SNPs"]):c += [x for x in hit["candidate_SNPs"].split(";") if x]
    out=[]
    for x in map(str.strip,c):
        if x and x not in out:out.append(x)
    return out

hits_by_tissue=defaultdict(list)
for h in hits:hits_by_tissue[canonical_tissue(h.get(tissue_col,"").strip())].append(h)
total_plotted=total_skipped=total_boxplots=0

for tissue,tissue_hits in hits_by_tissue.items():
    print(f"\n{'='*60}\nTISSUE: {tissue}\nPASSING HITS: {len(tissue_hits)}\n{'='*60}")
    qtl=ROOT/tissue/TYPE
    if TYPE in {"eQTL","circ_eQTL"}:pheno_file,loc_file=qtl/f"expression_{tissue}.txt",qtl/"gene_location.txt"
    elif TYPE=="sQTL":pheno_file,loc_file=qtl/f"splicing_{tissue}.txt",qtl/"splicing_location.txt"
    else:pheno_file,loc_file=qtl/f"circ_{tissue}.txt",qtl/"circ_location.txt"
    snp_file,snp_loc_file=qtl/f"snp_{tissue}.txt",qtl/"snp_location.txt";cov_file=qtl/f"covariates_{tissue}_encoded.txt"
    required=[pheno_file,loc_file,snp_file,snp_loc_file]
    if TYPE=="circ_eQTL":
        rpm_file=qtl/f"circ_rpm_{tissue}.txt";selected_file=qtl/"selected_samples.tsv"
        required += [rpm_file,selected_file]
    miss=[p for p in required if not p.exists()]
    if miss:print("SKIP TISSUE: missing "+", ".join(map(str,miss)));total_skipped+=len(tissue_hits);continue
    processed_roots_for_tissue(ROOT,tissue);status=load_covariate_status(cov_file)
    events={h[event_col].strip() for h in tissue_hits if not missing(h.get(event_col))};snps={s for h in tissue_hits for s in possible_snps(h)}
    pheno_rows=extract_matrix_rows(pheno_file,events);geno_rows=extract_matrix_rows(snp_file,snps)
    if TYPE=="circ_eQTL":
        rpm_rows=extract_matrix_rows(rpm_file,events)
        with open(selected_file,newline="") as f:
            selected_reader=csv.DictReader(f,delimiter="\t")
            if not {"donor_id","sample_id"}.issubset(selected_reader.fieldnames or []):
                sys.exit(f"ERROR: {selected_file} requires donor_id and sample_id")
            selected={}
            for row in selected_reader:
                donor,sample=row["donor_id"].strip(),row["sample_id"].strip()
                if donor in selected and norm(selected[donor])!=norm(sample):
                    sys.exit(f"ERROR: multiple selected samples for {donor} in {selected_file}")
                selected[donor]=sample
    trait_locs=load_location_table(loc_file,"trait",TYPE);snp_locs=load_location_table(snp_loc_file,"snp",TYPE)

    for i,hit in enumerate(tissue_hits,1):
        event=hit[event_col].strip();pp_h4=hit.get("max_PP_H4","");symbol=hit.get(symbol_col,"").strip() if symbol_col else ""
        print(f"\n[{i}/{len(tissue_hits)}] {tissue} | {event} | {symbol}")
        phenotype=find_matrix_row(pheno_rows,event)
        if phenotype is None:print("  SKIP: phenotype missing");total_skipped+=1;continue
        raw_rpm=None
        if TYPE=="circ_eQTL":
            raw_rpm=rpm_rows.get(event)
            if raw_rpm is None:print("  SKIP: exact circ ID absent from raw RPM matrix");total_skipped+=1;continue
            bad=[qid for qid in phenotype if qid not in selected or not selected[qid] or as_float(raw_rpm.get(qid)) is None or float(raw_rpm[qid])<0]
            if bad:print(f"  SKIP: missing selected sample or invalid raw RPM for {len(bad)} phenotype donors");total_skipped+=1;continue
        snp=genotype=None
        for candidate in possible_snps(hit):
            if candidate in geno_rows:snp,genotype=candidate,geno_rows[candidate];break
        if genotype is None or snp not in snp_locs:print("  SKIP: genotype/SNP location missing");total_skipped+=1;continue
        snp_chr,snp_pos=snp_locs[snp];allele=raw_alleles.get((chr_key(snp_chr),int(snp_pos)))
        if allele is None or allele["counted"] not in {allele["ref"],allele["alt"]}:print("  SKIP: allele orientation unavailable");total_skipped+=1;continue
        flip=allele["counted"]==allele["ref"];loc=parse_event_coordinates(event) if TYPE=="circ_eQTL" else find_location(trait_locs,event)
        if loc is None:print("  SKIP: coordinates unavailable");total_skipped+=1;continue
        print(f"  SNP alleles: REF={allele['ref']} ALT={allele['alt']} PLINK_counted={allele['counted']} | {'flipping to ALT dosage' if flip else 'already ALT dosage'}")
        chrom=loc["chr"];tstart,tend=sorted([loc["start"],loc["end"]]);strand=loc.get("strand")
        flank=min(int(max(FLANK_MIN,math.ceil(max(1,tend-tstart)*FLANK_FRAC))),FLANK_MAX);pstart,pend=max(0,tstart-flank),tend+flank
        nb=int(min(MAX_TRACK_POINTS,max(50,math.ceil((pend-pstart)/BIGWIG_BIN_BP))));binbp=(pend-pstart)/nb
        x=np.linspace(pstart,pend,nb,endpoint=False)+binbp/2;exons=load_reference_exons(GTF,chrom,pstart,pend)
        n_transcripts=len({e["transcript_id"] for e in exons if e["transcript_id"]})
        print(f"  GENCODE transcripts in plotted region: {n_transcripts}")
        label=symbol if symbol else event;base=safe_name(f"{tissue}__{TYPE}__{label}__{snp}")
        if TYPE=="circ_eQTL":
            # Full BSJ identity (including strand) prevents same-host filename collisions.
            base=safe_name(f"{tissue}__{TYPE}__{event}__{snp}")[:150]+"__"+hashlib.sha256(f"{tissue}|{event}|{snp}".encode()).hexdigest()[:16]
        boxout=OUTDIR/f"{base}.boxplot.pdf";pdf=OUTDIR/f"{base}.density.pdf";png=OUTDIR/f"{base}.density.png";svg=OUTDIR/f"{base}.density.svg";tsv=OUTDIR/f"{base}.subjects.tsv"

        records=[]
        for qid in phenotype.keys()&genotype.keys():
            g=genotype_group(genotype[qid],flip)
            try:ph=float(phenotype[qid])
            except:continue
            if g is None or not np.isfinite(ph):continue
            if TYPE=="circ_eQTL":
                # Use precisely the RNA sample selected during phenotype preparation.
                sid=selected[qid];subj=qid
            else:
                meta=metadata_for_qtl_id(qid,tissue)
                if not meta:continue
                sid=str(meta.get(sample_col,"")).strip();subj=str(meta.get(subject_col,"")).strip()
            sd=resolve_sample_dir(ROOT,tissue,sid) if sid else None
            if sd is None:
                print(f"    SAMPLE DIR MISSING: {qid} | {sid} | tissue={tissue}")
                continue
            records.append({"qtl_id":qid,"subject_id":subj,"sample_id":sid,"sample_fs_id":sd.name,"genotype":g,"genotype_raw":genotype[qid],"phenotype":ph,"is_als":status.get(qid,0),
                            "sample_dir":sd,"primary_reads":None,"plus_track":None,"minus_track":None,"plus_source":None,"minus_source":None,"circ_reads":None,"circ_rpm":None,"linear_reads":None,
                            "circ_percent_file":None,"circ_quant_source":None,"circ_match_type":None,"circ_matched_start":None,"circ_matched_end":None})
        print(f"  Candidate subjects before track QC: {len(records)}")
        for r in records:
            if TYPE=="circ_eQTL":
                r["circ_rpm"]=float(raw_rpm[r["qtl_id"]])
                r["circ_quant_source"]=str(rpm_file);r["circ_match_type"]="exact_matrix_ID"
                r["circ_matched_start"],r["circ_matched_end"]=tstart,tend
            if TYPE in {"sQTL","cQTL"}:r["primary_reads"]=primary_read_count(r["sample_dir"])
            r["plus_track"],r["minus_track"],r["plus_source"],r["minus_source"]=bigwig_strands(r["sample_dir"],chrom,pstart,pend,nb)
            if TYPE=="cQTL":
                qf=circ_quant_file(r["sample_dir"]);r["circ_quant_source"]=str(qf);q=read_circ_quant(qf,chrom,tstart,tend,strand,CIRC_COORD_TOL)
                if q:
                    r["circ_reads"],r["linear_reads"],r["circ_percent_file"]=q["circ_reads"],q["linear_reads"],q["circ_percent"]
                    r["circ_match_type"],r["circ_matched_start"],r["circ_matched_end"]=q["match_type"],q["matched_start"],q["matched_end"];r["circ_rpm"]=rpm_value(q["circ_reads"],r["primary_reads"])
        complete=[r for r in records if r["plus_track"] is not None and r["minus_track"] is not None and len(r["plus_track"])==nb and len(r["minus_track"])==nb and np.all(np.isfinite(r["plus_track"])) and np.all(np.isfinite(r["minus_track"]))]
        if not complete:print("  SKIP: no subjects with usable strand-specific tracks");total_skipped+=1;continue
        groups={g:[r for r in complete if r["genotype"]==g] for g in [0,1,2]};counts={g:len(groups[g]) for g in groups}
        print(f"  Matched plotted subjects: {len(complete)} | Ref/Ref={counts[0]} | Het={counts[1]} | Hom Alt={counts[2]}")
        if make_matched_boxplot(complete,boxout,event,snp,symbol,allele,TYPE):total_boxplots+=1
        plus={g:np.mean(np.vstack([r["plus_track"] for r in groups[g]]),0) if groups[g] else np.zeros(nb) for g in groups}
        minus={g:np.mean(np.vstack([r["minus_track"] for r in groups[g]]),0) if groups[g] else np.zeros(nb) for g in groups}
        junctions={}
        if TYPE in {"sQTL","cQTL"}:
            for g in groups:junctions[g],_=mean_junction_rpm(groups[g],chrom,pstart,pend,MIN_MEAN_JUNCTION_RPM,MAX_JUNCTIONS)
        circ={g:0. for g in groups}
        if TYPE=="circ_eQTL":
            for g in groups:
                if groups[g]:circ[g]=float(np.mean([r["circ_rpm"] for r in groups[g]]))
        if TYPE=="cQTL":
            for g in groups:
                valid=[r for r in groups[g] if r["primary_reads"]]
                if valid:circ[g]=sum(0 if r["circ_rpm"] is None else r["circ_rpm"] for r in valid)/len(valid)

        with open(tsv,"w",newline="") as f:
            w=csv.writer(f,delimiter="\t");w.writerow(["tissue","qtl_type","event","gene_symbol","snp","ref","alt","plink_counted_allele","qtl_id","subject_id","metadata_sample_id","filesystem_sample_id","genotype_group","genotype_label","genotype_raw","phenotype","is_als","primary_reads","plus_track_source","minus_track_source","circ_reads","circ_rpm","linear_reads","circ_percent_file","circ_match_type","circ_matched_start","circ_matched_end","circ_quant_source"])
            for r in complete:w.writerow([tissue,TYPE,event,symbol,snp,allele["ref"],allele["alt"],allele["counted"],r["qtl_id"],r["subject_id"],r["sample_id"],r["sample_fs_id"],r["genotype"],GENOTYPE_LABELS[r["genotype"]],r["genotype_raw"],r["phenotype"],r["is_als"],r["primary_reads"] or "",r["plus_source"],r["minus_source"],r["circ_reads"] or "","" if r["circ_rpm"] is None else r["circ_rpm"],r["linear_reads"] or "","" if r["circ_percent_file"] is None else r["circ_percent_file"],r["circ_match_type"] or "",r["circ_matched_start"] or "",r["circ_matched_end"] or "",r["circ_quant_source"] or ""])

        if TYPE!="cQTL" and not OVERWRITE and pdf.exists() and png.exists() and svg.exists():print("  DENSITY EXISTS; boxplot/TSV refreshed, density skipped");continue
        cmax=max(max(float(np.max(plus[g])),float(np.max(minus[g]))) for g in groups if groups[g]);cmax=1 if cmax<=0 else cmax*1.05
        linear=TYPE in {"sQTL","cQTL"} and any(junctions.get(g) for g in groups);circular=TYPE in {"cQTL","circ_eQTL"} and any(circ[g]>0 for g in groups)
        extent=cmax*(2.05 if circular else 1.80 if linear else 1.08);width=binbp*.92
        ref_ratio=max(.72,min(4.0,.18*max(1,n_transcripts)))
        fig_height=max(10,9+.22*n_transcripts)
        fig=plt.figure(figsize=(16,fig_height));gs=fig.add_gridspec(4,1,height_ratios=[1,1,1,ref_ratio],hspace=.10);axes=[]
        for g in [0,1,2]:
            ax=fig.add_subplot(gs[g,0]);axes.append(ax)
            if groups[g]:
                ax.bar(x,plus[g],width=width,color=TRACK_COLORS[g],alpha=.82,linewidth=0);ax.bar(x,-minus[g],width=width,color=TRACK_COLORS[g],alpha=.55,linewidth=0)
                if TYPE in {"sQTL","cQTL"}:draw_linear_arcs(ax,junctions[g],x,plus[g],minus[g],cmax,TRACK_COLORS[g])
                if TYPE in {"cQTL","circ_eQTL"}:draw_circ_arc(ax,tstart,tend,circ[g],cmax,TRACK_COLORS[g],strand)
            else:ax.text(.5,.5,"No subjects",transform=ax.transAxes,ha="center",va="center",fontsize=12)
            ax.axhline(0,color="#555555",lw=.8);ax.axvspan(tstart,tend,alpha=.10,color="#999999")
            if pstart<=snp_pos<=pend:ax.axvline(snp_pos,ls="--",lw=1.3,color="black",alpha=.8)
            ax.set_xlim(pstart,pend);ax.set_ylim(-extent,extent);ax.yaxis.set_major_formatter(FuncFormatter(density_formatter))
            ax.set_title(genotype_panel_title(groups[g],g,TYPE),loc="left",fontsize=13,color=TRACK_COLORS[g],y=.84,pad=0)
            ax.tick_params(axis="x",labelbottom=False);ax.tick_params(axis="y",labelsize=11);ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
            ax.text(.995,.94,"Forward (+)",transform=ax.transAxes,ha="right",va="top",fontsize=10);ax.text(.995,.06,"Reverse (-)",transform=ax.transAxes,ha="right",va="bottom",fontsize=10)
        ref=fig.add_subplot(gs[3,0],sharex=axes[0]);draw_reference_track(ref,exons,pstart,pend);ref.axvspan(tstart,tend,alpha=.10,color="#999999")
        if pstart<=snp_pos<=pend:ref.axvline(snp_pos,ls="--",lw=1,color="black",alpha=.7)
        ref.set_xlabel(f"{chrom} genomic position",fontsize=14,labelpad=10)
        title_event=f"{symbol} ({event})" if symbol and symbol!=event else event
        title=f"{tissue_display(tissue)} | {TYPE} | {title_event} | {snp} | REF={allele['ref']} ALT={allele['alt']}"
        if pp_h4 and not missing(pp_h4):
            try:title+=f" | PP.H4={float(pp_h4):.4f}"
            except:title+=f" | PP.H4={pp_h4}"
        fig.suptitle(title,fontsize=17,y=.995);fig.text(.008,.57,"Mean RPM-normalized RNA-seq read density",rotation=90,va="center",fontsize=15)
        footer=(f"Trait: {chrom}:{tstart:,}-{tend:,} | flank: {flank:,} bp each side | shaded region = trait locus | dashed vertical line = {snp} | coverage bin width ≈ {round(binbp)} bp\n"
                "Coverage uses the same RPM normalization as the source BigWigs (1e6 / primary alignments); + strand above 0, - strand below 0")
        if TYPE=="cQTL":footer+="\nPanel mean circ = 100 × mean per-subject circ_percent; missing/invalid quantification excluded (quant N shown when different).\nLinear splice junctions = solid arcs | circular/back-splice junction = dashed arc | labels = mean junction RPM per subject"
        elif TYPE=="circ_eQTL":footer+="\nDashed arc and panel label = mean BSJ RPM from the prepared circ_rpm matrix, including zeros; selected RNA samples only.\nBoxplots use the rank inverse-normal phenotype; coverage retains source BigWig normalization. Linear junction arcs are not shown."
        elif TYPE=="sQTL":footer+="\nArc labels = mean LeafCutter junction RPM per subject"
        fig.text(.5,.004,footer,ha="center",va="bottom",fontsize=9.5)
        fig.savefig(pdf,bbox_inches="tight");fig.savefig(png,dpi=DPI,bbox_inches="tight");fig.savefig(svg,bbox_inches="tight");plt.close(fig)
        print(f"  BOXPLOT PDF : {boxout}\n  DENSITY PDF : {pdf}\n  DENSITY PNG : {png}\n  DENSITY SVG : {svg}\n  SUBJECT TSV : {tsv}");total_plotted+=1

print(f"\n{'='*60}\nDONE\nQTL type             : {TYPE}\nDensity plots        : {total_plotted}\nBoxplots regenerated : {total_boxplots}\nSkipped hits         : {total_skipped}\nOutput               : {OUTDIR}\n{'='*60}")

PY
