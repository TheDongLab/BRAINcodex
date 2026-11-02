#!/usr/bin/env bash
#SBATCH --job-name=sQTL_SMR_HEIDI
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/MR/sQTL_SMR_HEIDI_5e8.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/MR/sQTL_SMR_HEIDI_5e8.err
#SBATCH --time=04:00:00
#SBATCH --mem=56G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
module load PLINK/1.9b_7.11-x86_64

QTL_LABEL="sQTL"
SMR="$HOME/donglab/pipelines/modules/smr/smr-1.4.2-linux-x86_64/smr"
EXPERIMENT="${EXPERIMENT:-target_ALS}"
[[ "$EXPERIMENT" == "target_ALS" || "$EXPERIMENT" == "answer_ALS" ]] || { echo "ERROR: EXPERIMENT must be target_ALS or answer_ALS" >&2; exit 1; }
BASE="$HOME/donglab/data/$EXPERIMENT"; GLOBAL="$BASE/MR/sQTL_SMR_HEIDI"
BFILE="$BASE/QTL/plink/joint_all_chrs_filtered_bed"; RAW="$BASE/QTL/plink/joint_all_chrs_matrixEQTL.raw"; BIM="$BFILE.bim"
GWAS_DIR="$HOME/donglab/references/GWAS/ALS"; GWAS="$GWAS_DIR/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"
GWAS_ORIG="$GWAS_DIR/GCST90027163_buildGRCh37.tsv.gz"; GWAS_MA="$GWAS_DIR/ALS_GRCh38_SMR.ma"
if [[ "$EXPERIMENT" == "answer_ALS" ]]; then
  TISSUES=(iPSC_Motor_Neuron)
else
  TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)
fi
[[ -z "${QTL_TISSUES:-}" ]] || read -r -a TISSUES <<< "$QTL_TISSUES"

PEQTL_SMR=5e-8; PEQTL_HEIDI=1.57e-3; HEIDI_PASS=0.05
LD_UPPER=0.90; LD_LOWER=0.05; HEIDI_MIN=3; HEIDI_MAX=20; CIS_WIND=2000

mkdir -p "$GLOBAL"
for f in "$SMR" "$BFILE.bed" "$BFILE.bim" "$BFILE.fam" "$RAW" "$GWAS"; do [[ -s "$f" ]] || { echo "ERROR: missing/empty: $f"; exit 1; }; done
command -v python3 >/dev/null || { echo "ERROR: python3 unavailable"; exit 1; }
python3 -c 'import pandas,numpy' >/dev/null 2>&1 || { echo "ERROR: pandas/numpy unavailable"; exit 1; }

export QTL_LABEL SMR BASE GLOBAL BFILE RAW BIM GWAS GWAS_ORIG GWAS_MA PEQTL_SMR PEQTL_HEIDI HEIDI_PASS
export LD_UPPER LD_LOWER HEIDI_MIN HEIDI_MAX CIS_WIND TISSUE_LIST="${TISSUES[*]}"

# ---------- Persistent GRCh38 GWAS .ma ----------
GWAS_VALID=0
if [[ -s "$GWAS_MA" ]]; then
  if python3 - "$GWAS_MA" <<'PY'
import sys,pandas as pd
f=sys.argv[1]
try:
    d=pd.read_csv(f,sep="\t",nrows=10000); req=["SNP","A1","A2","freq","b","se","p","n"]
    assert list(d.columns)==req and len(d)>0
    assert d.SNP.notna().all() and d.A1.notna().all() and d.A2.notna().all()
    se=pd.to_numeric(d.se,errors="coerce"); p=pd.to_numeric(d.p,errors="coerce")
    assert se.notna().all() and (se>0).all() and p.notna().all() and ((p>=0)&(p<=1)).all()
except Exception as e:
    print("INVALID GWAS .ma:",e); sys.exit(1)
PY
  then GWAS_VALID=1; fi
fi

if [[ "$GWAS_VALID" -eq 1 ]]; then
  echo "Reusing validated GWAS SMR file: $GWAS_MA"
else
  [[ -s "$GWAS_ORIG" ]] || { echo "ERROR: cannot rebuild GWAS .ma; missing $GWAS_ORIG"; exit 1; }
  echo "Creating GRCh38 GWAS SMR file: $GWAS_MA"; rm -f "$GWAS_MA"
  python3 <<'PY'
import os,pandas as pd,numpy as np
g=os.environ["GWAS"]; go=os.environ["GWAS_ORIG"]; out=os.environ["GWAS_MA"]
want=["hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value","hm_effect_allele_frequency"]
d=pd.read_csv(g,sep="\t",compression="gzip",usecols=lambda x:x in want,dtype=str)
need=["hm_rsid","hm_effect_allele","hm_other_allele","hm_beta","standard_error","p_value"]; miss=set(need)-set(d.columns)
if miss: raise RuntimeError(f"GWAS missing columns: {sorted(miss)}")
for c in ["hm_beta","standard_error","p_value"]: d[c]=pd.to_numeric(d[c],errors="coerce")
d["freq"]=pd.to_numeric(d["hm_effect_allele_frequency"],errors="coerce") if "hm_effect_allele_frequency" in d else np.nan
n=pd.read_csv(go,sep="\t",compression="gzip",usecols=["N_effective"]); N=int(round(pd.to_numeric(n.N_effective,errors="coerce").median()))
d=d.rename(columns={"hm_rsid":"SNP","hm_effect_allele":"A1","hm_other_allele":"A2","hm_beta":"b","standard_error":"se","p_value":"p"})
d=d[d.SNP.notna()&d.A1.notna()&d.A2.notna()&d.b.notna()&d.se.notna()&(d.se>0)&d.p.notna()].copy()
d["A1"]=d.A1.str.upper(); d["A2"]=d.A2.str.upper(); d=d[d.A1.isin(list("ACGT"))&d.A2.isin(list("ACGT"))].drop_duplicates("SNP"); d["n"]=N
d[["SNP","A1","A2","freq","b","se","p","n"]].to_csv(out,sep="\t",index=False,na_rep="NA")
print(f"Created {out}: {len(d):,} SNPs; N={N:,}")
PY
fi

# ---------- Per tissue ----------
for TISSUE in "${TISSUES[@]}"; do
  echo; echo "======================================================================"; echo "$TISSUE"; echo "======================================================================"

  SQTL="$BASE/$TISSUE/sQTL"; RESULTS="$SQTL/results"; OUT="$BASE/$TISSUE/MR/sQTL_SMR_HEIDI"
  INPUT="$OUT/input"; BESDDIR="$OUT/besd"; RESDIR="$OUT/results"; LDDIR="$OUT/ld_reference"
  FULL="$RESULTS/${TISSUE}_sQTL.full_annotated.txt"; SPLICELOC="$SQTL/splicing_location.txt"
  SNPLOC="$SQTL/snp_location.txt"; COV="$SQTL/covariates_${TISSUE}_encoded.txt"
  MATEQTL="$INPUT/${TISSUE}_MatrixEQTL_junction.txt"; KEEP="$INPUT/plink_keep.txt"; FREQ="$INPUT/${TISSUE}_freq"
  ESI_POOL="$INPUT/${TISSUE}_update.esi"; EPI_POOL="$INPUT/${TISSUE}_update.epi"; FREQ_POOL="$INPUT/${TISSUE}_update.freq"
  BESD="$BESDDIR/${TISSUE}_sQTL"; SMROUT="$RESDIR/${TISSUE}_sQTL_SMR_HEIDI"
  ESI_OK="$BESD.esi_native_updated.ok"; EPI_OK="$BESD.epi_native_updated.ok"; FREQ_OK="$BESD.freq_native_updated.ok"
  LD_BFILE="$LDDIR/${TISSUE}_sQTL_SMR_LD"; LD_MAP="$LDDIR/${TISSUE}_old_to_rsid.txt"; LD_EXTRACT="$LDDIR/${TISSUE}_extract.txt"; LD_OK="$LDDIR/${TISSUE}_sQTL_SMR_LD.ok"

  mkdir -p "$INPUT" "$BESDDIR" "$RESDIR" "$LDDIR"
  for f in "$FULL" "$SPLICELOC" "$SNPLOC" "$COV"; do [[ -s "$f" ]] || { echo "ERROR: missing/empty: $f"; exit 1; }; done
  export TISSUE SQTL RESULTS OUT INPUT BESDDIR RESDIR LDDIR FULL SPLICELOC SNPLOC COV MATEQTL KEEP FREQ
  export ESI_POOL EPI_POOL FREQ_POOL BESD SMROUT ESI_OK EPI_OK FREQ_OK LD_BFILE LD_MAP LD_EXTRACT LD_OK

  # ---------- Tissue subjects ----------
  if [[ ! -s "$KEEP" ]]; then
    python3 <<'PY'
import os,pandas as pd
fam=pd.read_csv(os.environ["BFILE"]+".fam",sep=r"\s+",header=None,usecols=[0,1],names=["FID","IID"],dtype=str)
cov=pd.read_csv(os.environ["COV"],sep="\t",nrows=0); subjects=set(map(str,list(cov.columns)[1:])); keep=fam[fam.IID.isin(subjects)]
if keep.empty: raise RuntimeError("No covariate subjects matched PLINK FAM")
keep.to_csv(os.environ["KEEP"],sep="\t",header=False,index=False)
print(f"{os.environ['TISSUE']}: sQTL/LD N={len(keep):,}")
PY
  fi
  N_SQTL=$(wc -l < "$KEEP"); echo "$TISSUE: sQTL/LD N=$N_SQTL"

  # ---------- MatrixEQTL -> BESD ----------
  if [[ -s "$BESD.besd" && -s "$BESD.esi" && -s "$BESD.epi" ]]; then
    echo "Reusing existing BESD: $BESD"
  else
    echo "Building BESD for $TISSUE"; rm -f "$BESD.besd" "$BESD.esi" "$BESD.epi" "$ESI_OK" "$EPI_OK" "$FREQ_OK" "$MATEQTL" "$MATEQTL.gz"

    python3 <<'PY'
import os,pandas as pd
T=os.environ["TISSUE"]; fullf=os.environ["FULL"]; locf=os.environ["SNPLOC"]; rawf=os.environ["RAW"]; bimf=os.environ["BIM"]; outf=os.environ["MATEQTL"]

with open(rawf) as fh: hdr=fh.readline().split()
meta={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}
rawids=set(c.rsplit("_",1)[0] for c in hdr if c not in meta and "_" in c)

bim=pd.read_csv(bimf,sep=r"\s+",header=None,names=["chr","rawid","cm","pos","a1","a2"],dtype=str)
bim=bim[bim.rawid.isin(rawids)].copy()
bim["key"]=bim.chr.str.replace("^chr","",regex=True).str.upper()+":"+bim.pos

loc=pd.read_csv(locf,sep="\t",dtype=str)
loc["key"]=loc.chr.str.replace("^chr","",regex=True).str.upper()+":"+loc.pos
available=set(loc[["snpid","key"]].drop_duplicates().merge(bim[["key"]].drop_duplicates(),on="key").snpid.astype(str))

use=["junction_id","snpid","beta","t-stat","p-value","FDR"]
first=True; total=written=duplicates=0

for ch in pd.read_csv(fullf,sep="\t",usecols=use,dtype=str,chunksize=250000):
    total+=len(ch); ch=ch[ch.snpid.isin(available)].copy()
    if ch.empty: continue

    for c in ["beta","t-stat","p-value","FDR"]: ch[c]=pd.to_numeric(ch[c],errors="coerce")
    ok=ch.beta.notna()&ch["t-stat"].notna()&(ch["t-stat"]!=0)&ch["p-value"].notna()&(ch["p-value"]>=0)&(ch["p-value"]<=1)&ch.FDR.notna()
    ch=ch[ok]

    out=ch[["snpid","junction_id","beta","t-stat","p-value","FDR"]].rename(columns={"snpid":"SNP","junction_id":"gene"})

    dup=out.duplicated(["SNP","gene"],keep=False)
    if dup.any():
        z=out.loc[dup,["SNP","gene"]].drop_duplicates()
        print(f"ERROR: {T}: duplicate SNP×junction mappings detected:")
        print(z.head(20).to_string(index=False)); duplicates+=len(z)

    out.to_csv(outf,sep="\t",index=False,mode="w" if first else "a",header=first)
    first=False; written+=len(out)

if duplicates: raise RuntimeError(f"{T}: {duplicates} duplicate SNP×junction mappings remain")
if not written: raise RuntimeError(f"{T}: no usable sQTL rows written for BESD")
print(f"{T}: source={total:,}; BESD input={written:,}")
PY

    "$SMR" --eqtl-summary "$MATEQTL" --matrix-eqtl-format --make-besd --add-n "$N_SQTL" --out "$BESD"
    [[ -s "$BESD.besd" && -s "$BESD.esi" && -s "$BESD.epi" ]] || { echo "ERROR: BESD creation failed"; exit 1; }
  fi

  # ---------- Native SMR metadata ----------
  if [[ ! -e "$ESI_OK" && -s "$BESD.esi.preupdate" ]]; then
    echo "Restoring pre-timeout ESI backup before native SMR update"
    cp -f "$BESD.esi.preupdate" "$BESD.esi"; rm -f "$BESD.esi.preupdate"
  fi

  if [[ ! -s "$FREQ.frq" ]]; then
    echo "Calculating tissue-specific allele frequencies"
    plink --bfile "$BFILE" --keep "$KEEP" --freq --threads "${SLURM_CPUS_PER_TASK:-1}" --out "$FREQ" >/dev/null
  fi

  if [[ ! -e "$ESI_OK" || ! -e "$EPI_OK" || ! -e "$FREQ_OK" ]]; then
    echo "Preparing native SMR metadata update files"
    python3 <<'PY'
import os,pandas as pd,numpy as np
T=os.environ["TISSUE"]; besd=os.environ["BESD"]; rawf=os.environ["RAW"]; bimf=os.environ["BIM"]; locf=os.environ["SNPLOC"]
traitf=os.environ["SPLICELOC"]; frqf=os.environ["FREQ"]+".frq"
esiout=os.environ["ESI_POOL"]; epiout=os.environ["EPI_POOL"]; freqout=os.environ["FREQ_POOL"]

esi=pd.read_csv(besd+".esi",sep=r"\s+",header=None,dtype=str)
epi=pd.read_csv(besd+".epi",sep=r"\s+",header=None,dtype=str)
besd_snps=pd.DataFrame({"snpid":esi.iloc[:,1].astype(str).unique()})
besd_probes=pd.DataFrame({"junction_id":epi.iloc[:,1].astype(str).unique()})

with open(rawf) as fh: hdr=fh.readline().split()
meta={"FID","IID","PAT","MAT","SEX","PHENOTYPE"}; rr=[]
for c in hdr:
    if c in meta or "_" not in c: continue
    vid,a=c.rsplit("_",1); rr.append((vid,a.upper()))
raw=pd.DataFrame(rr,columns=["rawid","effect"]).drop_duplicates("rawid")

bim=pd.read_csv(bimf,sep=r"\s+",header=None,names=["chr_bim","rawid","cm","pos_bim","a1","a2"],dtype=str)
bim["a1"]=bim.a1.str.upper(); bim["a2"]=bim.a2.str.upper()
bim["key"]=bim.chr_bim.str.replace("^chr","",regex=True).str.upper()+":"+bim.pos_bim
bim=bim.merge(raw,on="rawid",how="inner")
bim["other"]=np.where(bim.effect==bim.a1,bim.a2,np.where(bim.effect==bim.a2,bim.a1,np.nan))

loc=pd.read_csv(locf,sep="\t",dtype=str)
loc["chr_clean"]=loc.chr.str.replace("^chr","",regex=True)
loc["key"]=loc.chr_clean.str.upper()+":"+loc.pos

m=besd_snps.merge(loc[["snpid","chr_clean","pos","key"]].drop_duplicates("snpid"),on="snpid",how="left")
m=m.merge(bim[["rawid","key","effect","other"]].drop_duplicates("key"),on="key",how="left")

if m[["chr_clean","pos","rawid","effect","other"]].isna().any().any():
    z=m[m[["chr_clean","pos","rawid","effect","other"]].isna().any(axis=1)]
    raise RuntimeError(f"{T}: {len(z):,}/{len(m):,} BESD SNPs missing metadata")
if m.snpid.duplicated().any(): raise RuntimeError(f"{T}: duplicate SNP IDs after metadata merge")

m[["chr_clean","snpid"]].assign(GD="0",BP=m["pos"],A1=m["effect"],A2=m["other"])[["chr_clean","snpid","GD","BP","A1","A2"]].to_csv(esiout,sep="\t",header=False,index=False)

frq=pd.read_csv(frqf,sep=r"\s+",dtype=str); frq["MAF"]=pd.to_numeric(frq.MAF,errors="coerce")
f=m[["snpid","rawid","effect","other"]].merge(frq[["SNP","A1","A2","MAF"]],left_on="rawid",right_on="SNP",how="left")
f["A1"]=f.A1.str.upper(); f["A2"]=f.A2.str.upper()
f["freq"]=np.where(f.effect==f.A1,f.MAF,np.where(f.effect==f.A2,1-f.MAF,np.nan))
if f.freq.isna().any(): raise RuntimeError(f"{T}: {f.freq.isna().sum():,} BESD SNPs lack usable allele frequency")
f[["snpid","effect","other","freq"]].to_csv(freqout,sep="\t",header=False,index=False)

p=pd.read_csv(traitf,sep="\t",dtype=str)
need={"geneid","chr","left","right"}
if not need.issubset(p.columns): raise RuntimeError(f"{T}: splicing_location.txt missing columns: {sorted(need-set(p.columns))}")

p=p.rename(columns={"geneid":"junction_id"})
p["chr_clean"]=p.chr.str.replace("^chr","",regex=True)
p["left"]=pd.to_numeric(p.left,errors="coerce"); p["right"]=pd.to_numeric(p.right,errors="coerce")
p["bp"]=((p.left+p.right)/2).round().astype("Int64")
p=p.drop_duplicates("junction_id")

p=besd_probes.merge(p[["junction_id","chr_clean","bp"]],on="junction_id",how="left")
p["symbol"]=p.junction_id; p["strand"]="NA"; p["GD"]="0"

if p[["chr_clean","bp"]].isna().any().any():
    z=p[p[["chr_clean","bp"]].isna().any(axis=1)]
    raise RuntimeError(f"{T}: {len(z):,}/{len(p):,} BESD junctions missing coordinates")
if p.junction_id.duplicated().any(): raise RuntimeError(f"{T}: duplicate junction probes after metadata merge")

p[["chr_clean","junction_id","GD","bp","symbol","strand"]].to_csv(epiout,sep="\t",header=False,index=False)
print(f"{T}: update pools ready: {len(m):,} SNPs; {len(p):,} junctions")
PY
  fi

  if [[ ! -e "$ESI_OK" ]]; then echo "Updating ESI with native SMR --update-esi"; "$SMR" --beqtl-summary "$BESD" --update-esi "$ESI_POOL"; touch "$ESI_OK"; else echo "ESI already updated; skipping"; fi
  if [[ ! -e "$EPI_OK" ]]; then echo "Updating EPI with native SMR --update-epi"; "$SMR" --beqtl-summary "$BESD" --update-epi "$EPI_POOL"; touch "$EPI_OK"; else echo "EPI already updated; skipping"; fi
  if [[ ! -e "$FREQ_OK" ]]; then echo "Updating effect-allele frequencies with native SMR --update-freq"; "$SMR" --beqtl-summary "$BESD" --update-freq "$FREQ_POOL"; touch "$FREQ_OK"; else echo "Frequencies already updated; skipping"; fi

  python3 <<'PY'
import os,pandas as pd
T=os.environ["TISSUE"]; b=os.environ["BESD"]
e=pd.read_csv(b+".esi",sep=r"\s+",header=None,dtype=str); p=pd.read_csv(b+".epi",sep=r"\s+",header=None,dtype=str)
if e.shape[1]<6 or p.shape[1]<6: raise RuntimeError(f"{T}: malformed ESI/EPI after update")
if e.iloc[:,1].duplicated().any() or p.iloc[:,1].duplicated().any(): raise RuntimeError(f"{T}: duplicate SNP/probe IDs after update")
if (e.iloc[:,0].isin(["0","NA","nan"])).any() or (e.iloc[:,3].isin(["0","NA","nan"])).any(): raise RuntimeError(f"{T}: incomplete ESI coordinates")
if (p.iloc[:,0].isin(["0","NA","nan"])).any() or (p.iloc[:,3].isin(["0","NA","nan"])).any(): raise RuntimeError(f"{T}: incomplete EPI coordinates")
print(f"{T}: metadata integrity OK: {len(e):,} SNPs; {len(p):,} junction probes")
PY

  rm -f "$ESI_POOL" "$EPI_POOL" "$FREQ_POOL"
  [[ -s "$MATEQTL" ]] && gzip -f "$MATEQTL"

  # ---------- rsID-matched tissue LD reference ----------
  if [[ -e "$LD_OK" && -s "$LD_BFILE.bed" && -s "$LD_BFILE.bim" && -s "$LD_BFILE.fam" ]]; then
    echo "Reusing rsID-matched SMR LD reference: $LD_BFILE"
  else
    echo "Building rsID-matched SMR LD reference"
    rm -f "$LD_OK" "$LD_BFILE".{bed,bim,fam,log,nosex} "$LDDIR/${TISSUE}_sQTL_SMR_pre".{bed,bim,fam,log,nosex}

    python3 <<'PY'
import os,pandas as pd
T=os.environ["TISSUE"]; bimf=os.environ["BIM"]; locf=os.environ["SNPLOC"]; besd=os.environ["BESD"]
mapout=os.environ["LD_MAP"]; extractout=os.environ["LD_EXTRACT"]

bim=pd.read_csv(bimf,sep=r"\s+",header=None,names=["chr","old_id","cm","pos","a1","a2"],dtype=str)
bim["key"]=bim.chr.str.replace("^chr","",regex=True).str.upper()+":"+bim.pos

loc=pd.read_csv(locf,sep="\t",dtype=str)
loc["key"]=loc.chr.str.replace("^chr","",regex=True).str.upper()+":"+loc.pos

esi=pd.read_csv(besd+".esi",sep=r"\s+",header=None,dtype=str)
besd_rs=set(esi.iloc[:,1].astype(str))

unique_old=set(bim.loc[~bim.old_id.duplicated(keep=False),"old_id"])
x=bim[bim.old_id.isin(unique_old)][["old_id","key"]].merge(loc[["snpid","key"]].drop_duplicates(),on="key",how="inner")
x=x[x.snpid.astype(str).isin(besd_rs)].drop_duplicates(["old_id","snpid"])

old_n=x.groupby("old_id").snpid.nunique(); rs_n=x.groupby("snpid").old_id.nunique()
x=x[x.old_id.isin(old_n[old_n==1].index)&x.snpid.isin(rs_n[rs_n==1].index)].drop_duplicates(["old_id","snpid"])

if x.empty: raise RuntimeError(f"{T}: no unambiguous PLINK ID -> rsID mappings")
if x.old_id.duplicated().any() or x.snpid.duplicated().any(): raise RuntimeError(f"{T}: non-unique LD rename map remained")

x[["old_id","snpid"]].to_csv(mapout,sep="\t",header=False,index=False)
x[["old_id"]].to_csv(extractout,sep="\t",header=False,index=False)
print(f"{T}: unambiguous PLINK -> rsID mappings: {len(x):,}")
PY

    PRE="$LDDIR/${TISSUE}_sQTL_SMR_pre"
    plink --bfile "$BFILE" --keep "$KEEP" --extract "$LD_EXTRACT" --make-bed --threads "${SLURM_CPUS_PER_TASK:-1}" --out "$PRE" >/dev/null
    plink --bfile "$PRE" --update-name "$LD_MAP" --make-bed --threads "${SLURM_CPUS_PER_TASK:-1}" --out "$LD_BFILE" >/dev/null
    rm -f "$PRE".{bed,bim,fam,log,nosex}

    python3 <<'PY'
import os,pandas as pd
T=os.environ["TISSUE"]; b=os.environ["LD_BFILE"]
d=pd.read_csv(b+".bim",sep=r"\s+",header=None,dtype=str); fam=pd.read_csv(b+".fam",sep=r"\s+",header=None,dtype=str)
if d.empty: raise RuntimeError(f"{T}: empty SMR LD reference")
if d.iloc[:,1].duplicated().any(): raise RuntimeError(f"{T}: duplicate rsIDs remain in SMR LD reference")
print(f"{T}: SMR LD reference ready: {len(d):,} unique rsID SNPs; N={len(fam):,}")
PY
    touch "$LD_OK"
  fi

  # ---------- Standard SMR + HEIDI ----------
  SINGLE_VALID=0
  if [[ -s "$SMROUT.smr" ]] && head -1 "$SMROUT.smr" | grep -q "p_SMR" && head -1 "$SMROUT.smr" | grep -q "p_HEIDI"; then SINGLE_VALID=1; fi

  if [[ "$SINGLE_VALID" -eq 1 ]]; then
    echo "Reusing single-SNP sQTL SMR + HEIDI result: $SMROUT.smr"
  else
    rm -f "$SMROUT.smr"; echo "Running single-SNP sQTL SMR + HEIDI for $TISSUE"
    "$SMR" --bfile "$LD_BFILE" --gwas-summary "$GWAS_MA" --beqtl-summary "$BESD" \
      --peqtl-smr "$PEQTL_SMR" --peqtl-heidi "$PEQTL_HEIDI" --ld-upper-limit "$LD_UPPER" --ld-lower-limit "$LD_LOWER" \
      --heidi-min-m "$HEIDI_MIN" --heidi-max-m "$HEIDI_MAX" --heidi-mtd 1 --cis-wind "$CIS_WIND" \
      --thread-num "${SLURM_CPUS_PER_TASK:-1}" --out "$SMROUT"
    [[ -s "$SMROUT.smr" ]] || { echo "ERROR: no single-SNP SMR output for $TISSUE"; exit 1; }
  fi

  # ---------- SMR-multi ----------
  MULTI_VALID=0
  if [[ -s "$SMROUT.msmr" ]] && head -1 "$SMROUT.msmr" | grep -q "p_SMR_multi"; then MULTI_VALID=1; fi

  if [[ "$MULTI_VALID" -eq 1 ]]; then
    echo "Reusing sQTL SMR-multi result: $SMROUT.msmr"
  else
    rm -f "$SMROUT.msmr" "$SMROUT.snps4msmr.list"; echo "Running sQTL SMR-multi for $TISSUE"
    "$SMR" --bfile "$LD_BFILE" --gwas-summary "$GWAS_MA" --beqtl-summary "$BESD" \
      --peqtl-smr "$PEQTL_SMR" --peqtl-heidi "$PEQTL_HEIDI" --ld-upper-limit "$LD_UPPER" --ld-lower-limit "$LD_LOWER" \
      --heidi-min-m "$HEIDI_MIN" --heidi-max-m "$HEIDI_MAX" --heidi-mtd 1 --cis-wind "$CIS_WIND" \
      --smr-multi --thread-num "${SLURM_CPUS_PER_TASK:-1}" --out "$SMROUT"
    [[ -s "$SMROUT.msmr" ]] || { echo "ERROR: no SMR-multi output for $TISSUE"; exit 1; }
  fi
done

# ---------- Combine .smr + .msmr; BH-FDR + HEIDI diagnostics ----------
python3 <<'PY'
import os,numpy as np,pandas as pd
BASE=os.environ["BASE"]; GLOBAL=os.environ["GLOBAL"]; Q=os.environ["QTL_LABEL"]
tissues=os.environ["TISSUE_LIST"].split(); HEIDI=float(os.environ["HEIDI_PASS"])

def bh(p):
    p=pd.to_numeric(p,errors="coerce"); out=pd.Series(np.nan,index=p.index,dtype=float); ok=p.notna()
    if not ok.any(): return out
    x=p[ok].clip(0,1); idx=x.sort_values().index; v=x.loc[idx].to_numpy(); m=len(v); a=np.empty(m); run=1.
    for i in range(m-1,-1,-1): run=min(run,v[i]*m/(i+1)); a[i]=min(run,1.)
    out.loc[idx]=a; return out

allx=[]; summary=[]
for t in tissues:
    prefix=f"{BASE}/{t}/MR/sQTL_SMR_HEIDI/results/{t}_sQTL_SMR_HEIDI"; sf=prefix+".smr"; mf=prefix+".msmr"
    if not os.path.exists(sf): raise RuntimeError(f"Missing single-SNP SMR result: {sf}")
    if not os.path.exists(mf): raise RuntimeError(f"Missing SMR-multi result: {mf}")

    s=pd.read_csv(sf,sep=r"\s+",dtype=str); m=pd.read_csv(mf,sep=r"\s+",dtype=str)
    if "probeID" not in s or "p_SMR" not in s or "p_HEIDI" not in s: raise RuntimeError(f"Malformed .smr: {sf}")
    if "probeID" not in m or "p_SMR_multi" not in m: raise RuntimeError(f"Malformed .msmr: {mf}")
    if s.probeID.duplicated().any(): raise RuntimeError(f"Duplicate probeID in {sf}")
    if m.probeID.duplicated().any(): raise RuntimeError(f"Duplicate probeID in {mf}")

    for c in ["p_SMR","p_HEIDI","b_SMR","se_SMR","p_GWAS","p_eQTL","nsnp_HEIDI"]:
        if c in s: s[c]=pd.to_numeric(s[c],errors="coerce")
    m["p_SMR_multi"]=pd.to_numeric(m["p_SMR_multi"],errors="coerce")

    d=s.merge(m[["probeID","p_SMR_multi"]],on="probeID",how="left",validate="one_to_one")
    if d.p_SMR_multi.isna().all(): raise RuntimeError(f"No p_SMR_multi values merged for {t}")

    if "p_eQTL" in d.columns: d=d.rename(columns={"p_eQTL":"p_sQTL"})

    d.insert(0,"tissue",t); d["FDR_SMR_tissue"]=bh(d.p_SMR); d["FDR_SMR_multi_tissue"]=bh(d.p_SMR_multi)
    d["HEIDI_tested"]=d.p_HEIDI.notna(); d["HEIDI_pass"]=d.HEIDI_tested&(d.p_HEIDI>=HEIDI)
    d["SMR_tissue_FDR_pass"]=d.FDR_SMR_tissue<0.05; d["SMR_multi_tissue_FDR_pass"]=d.FDR_SMR_multi_tissue<0.05
    allx.append(d)

    summary.append({
      "tissue":t,"SMR_tests":len(d),
      "SMR_nominal_P0.05":int((d.p_SMR<.05).sum()),"SMR_tissue_FDR0.05":int((d.FDR_SMR_tissue<.05).sum()),
      "SMR_multi_nominal_P0.05":int((d.p_SMR_multi<.05).sum()),"SMR_multi_tissue_FDR0.05":int((d.FDR_SMR_multi_tissue<.05).sum()),
      "HEIDI_tested":int(d.HEIDI_tested.sum()),"HEIDI_P_ge_0.05":int(d.HEIDI_pass.sum())
    })

x=pd.concat(allx,ignore_index=True)
x["FDR_SMR_global"]=bh(x.p_SMR); x["FDR_SMR_multi_global"]=bh(x.p_SMR_multi)
x["SMR_global_FDR_pass"]=x.FDR_SMR_global<0.05; x["SMR_multi_global_FDR_pass"]=x.FDR_SMR_multi_global<0.05
x["SMR_HEIDI_pass"]=x.SMR_global_FDR_pass&x.HEIDI_pass
x["SMR_MULTI_HEIDI_pass"]=x.SMR_multi_global_FDR_pass&x.HEIDI_pass

for t in tissues:
    out=f"{BASE}/{t}/MR/sQTL_SMR_HEIDI"
    d=x[x.tissue==t].sort_values(["FDR_SMR_multi_global","p_SMR_multi","p_SMR"])
    d.to_csv(f"{out}/{t}_{Q}_SMR_HEIDI_all.tsv.gz",sep="\t",index=False,compression="gzip")
    d[d.SMR_HEIDI_pass].to_csv(f"{out}/{t}_{Q}_SMR_HEIDI_pass.tsv",sep="\t",index=False)
    d[d.SMR_MULTI_HEIDI_pass].to_csv(f"{out}/{t}_{Q}_SMR_MULTI_HEIDI_pass.tsv",sep="\t",index=False)

s=pd.DataFrame(summary)
for i,r in s.iterrows():
    d=x[x.tissue==r.tissue]
    s.loc[i,"SMR_global_FDR0.05"]=int((d.FDR_SMR_global<.05).sum())
    s.loc[i,"SMR_global_FDR0.05_and_HEIDI_pass"]=int(d.SMR_HEIDI_pass.sum())
    s.loc[i,"SMR_multi_global_FDR0.05"]=int((d.FDR_SMR_multi_global<.05).sum())
    s.loc[i,"SMR_multi_global_FDR0.05_and_HEIDI_pass"]=int(d.SMR_MULTI_HEIDI_pass.sum())

s.to_csv(f"{GLOBAL}/{Q}_SMR_HEIDI_summary_all_tissues.tsv",sep="\t",index=False)
x.to_csv(f"{GLOBAL}/all_tissues_{Q}_SMR_HEIDI.tsv.gz",sep="\t",index=False,compression="gzip")
x[x.SMR_HEIDI_pass].to_csv(f"{GLOBAL}/all_tissues_{Q}_SMR_HEIDI_pass.tsv",sep="\t",index=False)
x[x.SMR_MULTI_HEIDI_pass].to_csv(f"{GLOBAL}/all_tissues_{Q}_SMR_MULTI_HEIDI_pass.tsv",sep="\t",index=False)

# ==============================================================
# REPORT ALL FDR-SIGNIFICANT sQTL SMR RESULTS
# ==============================================================

x["SMR_tissue_FDR_pass"]=x["FDR_SMR_tissue"]<0.05
x["SMR_multi_tissue_FDR_pass"]=x["FDR_SMR_multi_tissue"]<0.05
x["SMR_global_FDR_pass"]=x["FDR_SMR_global"]<0.05
x["SMR_multi_global_FDR_pass"]=x["FDR_SMR_multi_global"]<0.05

x["SMR_tissue_HEIDI_pass"]=x["SMR_tissue_FDR_pass"]&x["HEIDI_pass"]
x["SMR_multi_tissue_HEIDI_pass"]=x["SMR_multi_tissue_FDR_pass"]&x["HEIDI_pass"]
x["SMR_global_HEIDI_pass"]=x["SMR_global_FDR_pass"]&x["HEIDI_pass"]
x["SMR_multi_global_HEIDI_pass"]=x["SMR_multi_global_FDR_pass"]&x["HEIDI_pass"]

hits=x[
    x["SMR_tissue_FDR_pass"] |
    x["SMR_multi_tissue_FDR_pass"] |
    x["SMR_global_FDR_pass"] |
    x["SMR_multi_global_FDR_pass"]
].copy()

hits=hits.sort_values(
    ["FDR_SMR_global","FDR_SMR_multi_global","FDR_SMR_tissue","FDR_SMR_multi_tissue","p_SMR"]
)

hits.to_csv(
    f"{GLOBAL}/all_tissues_{Q}_SMR_any_FDR_pass_with_HEIDI.tsv",
    sep="\t",index=False
)

print("\n======================================================================")
print("sQTL SMR + SMR-MULTI + HEIDI SUMMARY")
print("======================================================================")
print(s.to_string(index=False))

print(f"\nTotal sQTL SMR tests: {len(x):,}")
print(f"Single-SNP sQTL SMR tissue FDR<0.05: {x.SMR_tissue_FDR_pass.sum():,}")
print(f"Single-SNP sQTL SMR tissue FDR<0.05 + HEIDI P>={HEIDI}: {x.SMR_tissue_HEIDI_pass.sum():,}")
print(f"sQTL SMR-multi tissue FDR<0.05: {x.SMR_multi_tissue_FDR_pass.sum():,}")
print(f"sQTL SMR-multi tissue FDR<0.05 + HEIDI P>={HEIDI}: {x.SMR_multi_tissue_HEIDI_pass.sum():,}")
print(f"Single-SNP sQTL SMR global FDR<0.05: {x.SMR_global_FDR_pass.sum():,}")
print(f"Single-SNP sQTL SMR global FDR<0.05 + HEIDI P>={HEIDI}: {x.SMR_global_HEIDI_pass.sum():,}")
print(f"sQTL SMR-multi global FDR<0.05: {x.SMR_multi_global_FDR_pass.sum():,}")
print(f"sQTL SMR-multi global FDR<0.05 + HEIDI P>={HEIDI}: {x.SMR_multi_global_HEIDI_pass.sum():,}")

print("\n======================================================================")
print("ALL sQTL SMR FDR<0.05 HITS")
print("======================================================================")

if len(hits):
    cols=[
        "tissue","probeID","Gene","topSNP",
        "p_sQTL","p_GWAS",
        "p_SMR","FDR_SMR_tissue","SMR_tissue_FDR_pass",
        "p_SMR_multi","FDR_SMR_multi_tissue","SMR_multi_tissue_FDR_pass",
        "FDR_SMR_global","SMR_global_FDR_pass",
        "FDR_SMR_multi_global","SMR_multi_global_FDR_pass",
        "p_HEIDI","nsnp_HEIDI","HEIDI_pass"
    ]
    print(hits[[c for c in cols if c in hits.columns]].to_string(index=False))
else:
    print("None")
PY

echo; echo "sQTL SMR + SMR-multi + HEIDI complete."; echo
column -t -s $'\t' "$GLOBAL/${QTL_LABEL}_SMR_HEIDI_summary_all_tissues.tsv"
