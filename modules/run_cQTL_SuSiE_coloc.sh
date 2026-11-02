#!/usr/bin/env bash
#SBATCH --job-name=cQTL_SuSiE_coloc
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/MR/cQTL_SuSiE_coloc.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/MR/cQTL_SuSiE_coloc.err
#SBATCH --time=08:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4

set -euo pipefail
module --force purge
module load PLINK/1.9b_7.11-x86_64
module load R

command -v plink >/dev/null || { echo "ERROR: PLINK unavailable"; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript unavailable"; exit 1; }
Rscript -e 'p<-c("data.table","susieR","coloc","ggplot2"); x<-p[!sapply(p,requireNamespace,quietly=TRUE)]; if(length(x)) stop("Missing R packages: ",paste(x,collapse=", "))'

EXPERIMENT="${EXPERIMENT:-target_ALS}"
[[ "$EXPERIMENT" == "target_ALS" || "$EXPERIMENT" == "answer_ALS" ]] || { echo "ERROR: EXPERIMENT must be target_ALS or answer_ALS" >&2; exit 1; }
BASE="$HOME/donglab/data/$EXPERIMENT"; BFILE="$BASE/QTL/plink/joint_all_chrs_filtered_bed"; RAW="$BASE/QTL/plink/joint_all_chrs_matrixEQTL.raw"
GWAS="$HOME/donglab/references/GWAS/ALS/harmonised/34873335-GCST90027163-MONDO_0004976.h.tsv.gz"; GWAS_ORIG="$HOME/donglab/references/GWAS/ALS/GCST90027163_buildGRCh37.tsv.gz"
if [[ "$EXPERIMENT" == "answer_ALS" ]]; then
  TISSUES=(iPSC_Motor_Neuron)
else
  TISSUES=(Cervical_Spinal_Cord Lumbar_Spinal_Cord Motor_Cortex Frontal_Cortex Cerebellum)
fi
[[ -z "${QTL_TISSUES:-}" ]] || read -r -a TISSUES <<< "$QTL_TISSUES"

CANDIDATE_MODE="global_fdr"   # global_fdr | tissue_fdr | nominal | all
H4_REFERENCE=0.80; MIN_SHARED_SNPS=10; MIN_SHARED_FRAC=0.50; SUSIE_L=10; MIN_F=10

export BASE BFILE RAW GWAS GWAS_ORIG CANDIDATE_MODE H4_REFERENCE MIN_SHARED_SNPS MIN_SHARED_FRAC SUSIE_L MIN_F TISSUE_LIST="${TISSUES[*]}"

Rscript --vanilla - <<'RS'
suppressPackageStartupMessages({library(data.table); library(susieR); library(coloc); library(ggplot2)})

BASE<-Sys.getenv("BASE"); BFILE<-Sys.getenv("BFILE"); RAW<-Sys.getenv("RAW"); GWAS<-Sys.getenv("GWAS"); GWAS_ORIG<-Sys.getenv("GWAS_ORIG")
TISSUES<-strsplit(Sys.getenv("TISSUE_LIST")," ")[[1]]; MODE<-Sys.getenv("CANDIDATE_MODE"); H4REF<-as.numeric(Sys.getenv("H4_REFERENCE"))
MINSHARED<-as.integer(Sys.getenv("MIN_SHARED_SNPS")); MINFRAC<-as.numeric(Sys.getenv("MIN_SHARED_FRAC")); L<-as.integer(Sys.getenv("SUSIE_L")); MINF<-as.numeric(Sys.getenv("MIN_F"))
comp<-c(A="T",T="A",C="G",G="C"); nch<-function(x) toupper(sub("^chr","",as.character(x),ignore.case=TRUE))

empty_cs<-function(file,status) fwrite(data.table(status=status,credible_set=NA_integer_,snps=NA_character_,n_snps=0L),file,sep="\t")
write_cs<-function(S,file,snps,label){
  z<-S$sets$cs
  if(is.null(z)||!length(z)){empty_cs(file,paste0("no_",label,"_credible_set")); return(0L)}
  fwrite(rbindlist(lapply(seq_along(z),function(i)data.table(status="OK",credible_set=i,snps=paste(snps[z[[i]]],collapse=";"),n_snps=length(z[[i]])))),file,sep="\t"); length(z)
}

# ==============================================================
# GENOTYPE ALLELE ORIENTATION
# ==============================================================
cat("Reading genotype allele orientation...\n")
bim<-fread(paste0(BFILE,".bim"),header=FALSE,col.names=c("chr","rawid","cm","pos","a1","a2")); bim[,key:=paste0(nch(chr),":",pos)]
hdr<-scan(RAW,what="",nlines=1,quiet=TRUE); hdr<-hdr[!hdr%in%c("FID","IID","PAT","MAT","SEX","PHENOTYPE")]
rawmap<-data.table(rawcol=hdr); rawmap[,`:=`(rawid=sub("_[^_]+$","",rawcol),effect_allele=toupper(sub("^.*_","",rawcol)))]
amap<-merge(rawmap,bim[,.(rawid,key,a1=toupper(a1),a2=toupper(a2))],by="rawid",all.x=TRUE)
amap[,other_allele:=fifelse(effect_allele==a1,a2,fifelse(effect_allele==a2,a1,NA_character_))]
amap<-unique(amap[!is.na(key),.(key,rawid,effect_allele,other_allele)],by="key"); rm(rawmap,hdr,bim); gc()

# ==============================================================
# GRCh38-HARMONIZED ALS GWAS
# ==============================================================
cat("Reading ALS GWAS...\n")
gwas<-fread(GWAS,select=c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"))
setnames(gwas,c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","standard_error","p_value"),
              c("snpid","gwas_chr","gwas_pos","gwas_oa","gwas_ea","gwas_beta","gwas_se","gwas_p"))
gwas[,`:=`(snpid=as.character(snpid),gwas_chr=nch(gwas_chr),gwas_pos=as.integer(gwas_pos),gwas_ea=toupper(gwas_ea),gwas_oa=toupper(gwas_oa),
           gwas_beta=as.numeric(gwas_beta),gwas_se=as.numeric(gwas_se),gwas_p=as.numeric(gwas_p))]
gwas<-unique(gwas[!is.na(snpid)],by="snpid")

gn<-fread(GWAS_ORIG,select=c("rsid","N_effective")); GWAS_N<-round(median(as.numeric(gn$N_effective),na.rm=TRUE)); rm(gn); gc()
fam<-fread(paste0(BFILE,".fam"),header=FALSE,select=1:2,col.names=c("FID","IID")); cat("GWAS effective N =",GWAS_N,"\n")
master<-list()

for(tissue in TISSUES){
  mrdir<-file.path(BASE,tissue,"MR"); snpfile<-file.path(mrdir,paste0(tissue,"_cQTL_SNP_results.tsv.gz"))
  outroot<-file.path(mrdir,"cQTL_SuSiE_coloc"); cqtldir<-file.path(BASE,tissue,"cQTL"); resultdir<-file.path(cqtldir,"results")
  fullfile<-file.path(resultdir,paste0(tissue,"_cQTL.full_annotated.txt")); circfile<-file.path(cqtldir,"circ_location.txt")
  locfile<-file.path(cqtldir,"snp_location.txt"); covfile<-file.path(cqtldir,paste0("covariates_",tissue,"_encoded.txt"))

  for(f in c(snpfile,fullfile,circfile,locfile,covfile)) if(!file.exists(f)||file.info(f)$size==0) stop("Missing/empty: ",f)
  if(dir.exists(outroot)) unlink(outroot,recursive=TRUE,force=TRUE); dir.create(outroot,recursive=TRUE,showWarnings=FALSE)

  # ============================================================
  # UNIQUE SNP CANDIDATES FROM CURRENT SNP-FIRST ANALYSIS
  # ============================================================
  snps<-fread(snpfile)
  for(c in intersect(c("gwas_p","FDR_ALS_tissue","FDR_ALS_global"),names(snps))) set(snps,j=c,value=as.numeric(snps[[c]]))

  cand_snps<-switch(MODE,
    global_fdr=snps[is.finite(FDR_ALS_global)&FDR_ALS_global<.05],
    tissue_fdr=snps[is.finite(FDR_ALS_tissue)&FDR_ALS_tissue<.05],
    nominal=snps[is.finite(gwas_p)&gwas_p<.05],
    all=copy(snps),
    stop("Unknown CANDIDATE_MODE: ",MODE))

  if(!nrow(cand_snps)){cat("\n",tissue,": zero candidate SNPs\n",sep=""); next}

  # ============================================================
  # READ DENSE cQTL TABLE
  # ============================================================
  full<-fread(fullfile)
  req<-c("circ_id","snpid","beta","t-stat","p-value","FDR","chr","pos","circ_chr","circ_start","circ_end")
  miss<-setdiff(req,names(full)); if(length(miss)) stop(tissue,": full_annotated missing: ",paste(miss,collapse=", "))

  full[,`:=`(
    circ_id=as.character(circ_id),snpid=as.character(snpid),
    beta_qtl=as.numeric(beta),t_qtl=as.numeric(`t-stat`),p_qtl=as.numeric(`p-value`),cqtl_FDR=as.numeric(FDR),
    circ_chr_key=nch(circ_chr),circ_start_key=as.integer(circ_start),circ_end_key=as.integer(circ_end)
  )]
  full[,F_statistic:=t_qtl^2]

  # ============================================================
  # ORIGINAL circRNA COORDINATE TABLE
  # ============================================================
  cloc<-fread(circfile)
  reqc<-c("circ_id","chr","left","right"); miss<-setdiff(reqc,names(cloc))
  if(length(miss)) stop(tissue,": circ_location.txt missing: ",paste(miss,collapse=", "))

  setnames(cloc,c("circ_id","chr","left","right"),c("circ_id","circ_chr","circ_start","circ_end"))
  cloc[,`:=`(
    circ_id=as.character(circ_id),circ_chr_key=nch(circ_chr),
    circ_start_key=as.integer(circ_start),circ_end_key=as.integer(circ_end)
  )]
  cloc<-unique(cloc[,.(circ_id,circ_chr_key,circ_start_key,circ_end_key)])

  if(cloc[,anyDuplicated(circ_id)]){
    dup<-cloc[,.(n_coords=.N),by=circ_id][n_coords>1]
    if(nrow(dup)) stop(tissue,": ",nrow(dup)," circRNA IDs have multiple coordinate records")
  }

  # ============================================================
  # CANDIDATE cQTL ASSOCIATIONS
  # ============================================================
  cand<-full[
    snpid%in%cand_snps$snpid &
    is.finite(cqtl_FDR)&cqtl_FDR<=.05 &
    is.finite(F_statistic)&F_statistic>=MINF,
    .(circ_id,snpid,beta_qtl,t_qtl,p_qtl,cqtl_FDR,F_statistic,
      circ_chr_key,circ_start_key,circ_end_key)
  ]

  cand<-unique(cand,by=c("circ_id","snpid"))
  cand<-merge(cand,cand_snps[,.(snpid,gwas_p,FDR_ALS_tissue,FDR_ALS_global)],by="snpid",all.x=TRUE)

  if(!nrow(cand)){cat("\n",tissue,": candidate SNPs had no qualifying cQTL relationships\n",sep=""); next}

  missing_circs<-setdiff(unique(cand$circ_id),cloc$circ_id)
  if(length(missing_circs)){
    fwrite(data.table(circ_id=missing_circs),file.path(outroot,"candidate_circRNAs_missing_location.tsv"),sep="\t")
    stop(tissue,": ",length(missing_circs)," candidate circRNAs missing from circ_location.txt")
  }

  circ_meta<-cand[,{
    z<-unique(.SD[,.(snpid,gwas_p,FDR_ALS_tissue,FDR_ALS_global)])[order(gwas_p)]
    .(
      n_candidate_SNPs=uniqueN(snpid),
      candidate_SNPs=paste(sort(unique(snpid)),collapse=";"),
      best_candidate_SNP=as.character(z$snpid[1]),
      best_ALS_p=as.numeric(z$gwas_p[1]),
      best_ALS_tissue_FDR=min(z$FDR_ALS_tissue,na.rm=TRUE),
      best_ALS_global_FDR=min(z$FDR_ALS_global,na.rm=TRUE)
    )
  },by=circ_id]

  circ_meta[!is.finite(best_ALS_tissue_FDR),best_ALS_tissue_FDR:=NA_real_]
  circ_meta[!is.finite(best_ALS_global_FDR),best_ALS_global_FDR:=NA_real_]
  circ_meta<-merge(circ_meta,cloc,by="circ_id",all.x=TRUE,validate="one-to-one")
  circs<-circ_meta$circ_id

  fwrite(cand,file.path(outroot,"candidate_SNP_circRNA_relationships.tsv"),sep="\t")

  cat("\n",strrep("=",70),"\n",tissue,": ",
      nrow(cand_snps)," unique candidate SNPs → ",
      nrow(cand)," SNP×circRNA relationships → ",
      length(circs)," candidate circRNA traits\n",
      strrep("=",70),"\n",sep="")

  # ============================================================
  # SNP LOCATIONS + TISSUE SUBJECTS
  # ============================================================
  loc<-fread(locfile)
  loc[,`:=`(snpid=as.character(snpid),qtl_chr=nch(chr),qtl_pos=as.integer(pos),key=paste0(nch(chr),":",pos))]
  loc<-unique(loc,by="snpid")

  covhdr<-names(fread(covfile,nrows=0)); subjects<-setdiff(covhdr,covhdr[1]); keep<-fam[IID%in%subjects]
  if(!nrow(keep)) stop(tissue,": no covariate subjects matched PLINK FAM")

  N_CQTL<-nrow(keep); keepfile<-file.path(outroot,"plink_keep.txt")
  fwrite(keep,keepfile,sep="\t",col.names=FALSE); cat("cQTL N =",N_CQTL,"\n")

  # ============================================================
  # SuSiE / COLOC PER circRNA
  # ============================================================
  for(circ in circs){
    cm<-circ_meta[circ_id==circ][1]; label<-circ
    safe<-gsub("[^A-Za-z0-9_.-]","_",circ)
    cat("  ",label," ... ",sep=""); cout<-file.path(outroot,safe); dir.create(cout,showWarnings=FALSE)

    q<-copy(full[circ_id==circ])

    if(!nrow(q)){
      cat("candidate circRNA missing from dense cQTL data\n")
      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status="candidate_circRNA_missing_from_current_cQTL",
        n_locus=0,n_shared=0,shared_fraction=0,n_candidate_SNPs=cm$n_candidate_SNPs,
        best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    q<-q[
      circ_chr_key==cm$circ_chr_key &
      circ_start_key==cm$circ_start_key &
      circ_end_key==cm$circ_end_key
    ]

    q[,`:=`(snpid=as.character(snpid),beta_qtl=as.numeric(beta),t_qtl=as.numeric(`t-stat`),p_qtl=as.numeric(`p-value`))]
    q[,se_qtl:=abs(beta_qtl/t_qtl)]
    q<-q[!is.na(snpid)&is.finite(beta_qtl)&is.finite(se_qtl)&se_qtl>0]
    q<-unique(q[,.(snpid,beta_qtl,se_qtl,p_qtl)],by="snpid")
    q<-merge(q,loc[,.(snpid,qtl_chr,qtl_pos,key)],by="snpid",all.x=TRUE)
    q<-merge(q,amap,by="key",all.x=TRUE)
    q<-unique(q,by="snpid"); n_locus<-nrow(q)

    if(n_locus==0){
      cat("no usable dense cQTL SNPs\n")
      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status="no_usable_dense_cQTL_snps",
        n_locus=0,n_shared=0,shared_fraction=0,n_candidate_SNPs=cm$n_candidate_SNPs,
        best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    if(n_locus<MINSHARED){
      cat("too few locus SNPs (",n_locus,")\n",sep="")
      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status="too_few_locus_snps",
        n_locus=n_locus,n_shared=0,shared_fraction=0,n_candidate_SNPs=cm$n_candidate_SNPs,
        best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    # ==========================================================
    # cQTL ↔ GRCh38 GWAS HARMONIZATION
    # ==========================================================
    x<-merge(q,gwas,by="snpid",all=FALSE)
    x<-x[qtl_chr==gwas_chr&qtl_pos==gwas_pos&!is.na(effect_allele)&!is.na(other_allele)]
    x[,pal:=paste0(effect_allele,other_allele)%in%c("AT","TA","CG","GC")]

    x[,relation:=fcase(
      effect_allele==gwas_ea&other_allele==gwas_oa,"exact",
      effect_allele==gwas_oa&other_allele==gwas_ea,"swap",
      effect_allele==comp[gwas_ea]&other_allele==comp[gwas_oa],"strand",
      effect_allele==comp[gwas_oa]&other_allele==comp[gwas_ea],"strand_swap",
      default="bad"
    )]

    x<-unique(x[!pal&relation!="bad"&is.finite(gwas_beta)&is.finite(gwas_se)&gwas_se>0],by="snpid")
    x[,gwas_beta_h:=fifelse(relation%in%c("swap","strand_swap"),-gwas_beta,gwas_beta)]
    shared_frac<-nrow(x)/n_locus

    if(nrow(x)<MINSHARED||shared_frac<MINFRAC){
      status<-if(nrow(x)<MINSHARED) "too_few_shared_snps" else "low_shared_coverage"
      cat(status,", shared=",nrow(x),"/",n_locus," (",round(100*shared_frac,1),"%)\n",sep="")

      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status,n_locus=n_locus,n_shared=nrow(x),
        shared_fraction=shared_frac,n_candidate_SNPs=cm$n_candidate_SNPs,
        best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    # ==========================================================
    # TISSUE-SPECIFIC LD
    # ==========================================================
    idsfile<-file.path(cout,"raw_ids.txt"); fwrite(x[,.(rawid)],idsfile,col.names=FALSE)
    pref<-file.path(cout,"ldgeno")

    system2("plink",c(
      "--bfile",BFILE,"--keep",keepfile,"--extract",idsfile,
      "--recode","A-transpose","--threads",Sys.getenv("SLURM_CPUS_PER_TASK","1"),
      "--out",pref
    ),stdout=TRUE,stderr=TRUE)

    trawfile<-paste0(pref,".traw")

    if(!file.exists(trawfile)){
      cat("PLINK failed\n")
      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status="PLINK_failed",
        n_locus=n_locus,n_shared=nrow(x),shared_fraction=shared_frac,
        n_candidate_SNPs=cm$n_candidate_SNPs,best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    tr<-fread(trawfile)
    geno_cols<-setdiff(names(tr),c("CHR","SNP","(C)M","POS","COUNTED","ALT"))

    if(nrow(tr)<MINSHARED||!length(geno_cols)){
      cat("too few PLINK SNPs\n")
      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status="too_few_PLINK_snps",
        n_locus=n_locus,n_shared=nrow(x),shared_fraction=shared_frac,
        n_candidate_SNPs=cm$n_candidate_SNPs,best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    G<-as.matrix(tr[,..geno_cols]); storage.mode(G)<-"numeric"
    orient<-x[match(tr$SNP,rawid),effect_allele]; flip<-toupper(tr$COUNTED)!=toupper(orient)
    G[flip,]<-2-G[flip,]
    G<-t(apply(G,1,function(v){v[is.na(v)]<-mean(v,na.rm=TRUE); v}))

    sdv<-apply(G,1,sd); ok<-is.finite(sdv)&sdv>0
    tr<-tr[ok]; G<-G[ok,,drop=FALSE]; x<-x[match(tr$SNP,rawid)]

    if(nrow(x)<MINSHARED){
      cat("too few variable SNPs\n")
      master[[length(master)+1]]<-data.table(
        tissue,circ_id=circ,status="too_few_variable_snps",
        n_locus=n_locus,n_shared=nrow(x),shared_fraction=nrow(x)/n_locus,
        n_candidate_SNPs=cm$n_candidate_SNPs,best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p)
      next
    }

    R<-cor(t(G))
    if(!is.matrix(R)||nrow(R)!=nrow(x)) stop(label,": malformed LD matrix")

    R[R>1]<-1; R[R< -1]<- -1; diag(R)<-1
    colnames(R)<-rownames(R)<-x$snpid

    # ==========================================================
    # SuSiE cQTL + SuSiE GWAS
    # ==========================================================
    d1<-list(beta=x$beta_qtl,varbeta=x$se_qtl^2,snp=x$snpid,position=x$qtl_pos,type="quant",sdY=1,N=N_CQTL,LD=R)
    d2<-list(beta=x$gwas_beta_h,varbeta=x$gwas_se^2,snp=x$snpid,position=x$qtl_pos,type="cc",N=GWAS_N,LD=R)

    invisible(coloc::check_dataset(d1,req="LD")); invisible(coloc::check_dataset(d2,req="LD"))

    z1<-d1$beta/sqrt(d1$varbeta); z2<-d2$beta/sqrt(d2$varbeta)
    s_cqtl<-tryCatch(susieR::estimate_s_rss(z1,R,n=N_CQTL),error=function(e) NA_real_)
    s_gwas<-tryCatch(susieR::estimate_s_rss(z2,R,n=GWAS_N),error=function(e) NA_real_)

    run_coloc_results<-function(S1,S2,x,cout){
      pip1<-if(is.null(S1$pip)) rep(NA_real_,nrow(x)) else S1$pip
      pip2<-if(is.null(S2$pip)) rep(NA_real_,nrow(x)) else S2$pip

      fwrite(data.table(snp=x$snpid,position=x$qtl_pos,cqtl_pip=pip1,gwas_pip=pip2),
             file.path(cout,"PIP.tsv"),sep="\t")

      ncs1<-write_cs(S1,file.path(cout,"cQTL_credible_sets.tsv"),x$snpid,"cQTL")
      ncs2<-write_cs(S2,file.path(cout,"GWAS_credible_sets.tsv"),x$snpid,"GWAS")

      top1<-if(any(is.finite(pip1))) which.max(pip1) else NA_integer_
      top2<-if(any(is.finite(pip2))) which.max(pip2) else NA_integer_

      ret<-function(status,h4=NA_real_) list(
        status=status,h4=h4,ncs1=ncs1,ncs2=ncs2,
        cqtl_top=if(!is.na(top1))x$snpid[top1] else NA_character_,
        cqtl_pip=if(!is.na(top1))pip1[top1] else NA_real_,
        gwas_top=if(!is.na(top2))x$snpid[top2] else NA_character_,
        gwas_pip=if(!is.na(top2))pip2[top2] else NA_real_
      )

      if(ncs1==0||ncs2==0){
        status<-if(ncs1==0&&ncs2==0) "no_cQTL_or_GWAS_credible_set" else if(ncs1==0) "no_cQTL_credible_set" else "no_GWAS_credible_set"
        fwrite(data.table(status=status),file.path(cout,"coloc_susie.tsv"),sep="\t")
        return(ret(status))
      }

      C<-coloc::coloc.susie(S1,S2)
      csum<-if(is.null(C$summary)) data.table() else as.data.table(C$summary)

      if(!nrow(csum)){
        fwrite(data.table(status="no_coloc_signal_pairs"),file.path(cout,"coloc_susie.tsv"),sep="\t")
        return(ret("no_coloc_signal_pairs"))
      }

      fwrite(csum,file.path(cout,"coloc_susie.tsv"),sep="\t")
      h4<-if("PP.H4.abf"%in%names(csum)) suppressWarnings(max(csum$PP.H4.abf,na.rm=TRUE)) else NA_real_
      if(!is.finite(h4)) h4<-NA_real_

      ret("OK",h4)
    }

    ans<-tryCatch({
      S1<-coloc::runsusie(d1,maxit=200,L=L)
      S2<-coloc::runsusie(d2,maxit=200,L=L)
      run_coloc_results(S1,S2,x,cout)
    },error=function(e) list(
      status=paste0("ERROR: ",conditionMessage(e)),h4=NA_real_,
      ncs1=NA_integer_,ncs2=NA_integer_,
      cqtl_top=NA_character_,cqtl_pip=NA_real_,
      gwas_top=NA_character_,gwas_pip=NA_real_
    ))

    fwrite(
      x[,.(snpid,qtl_chr,qtl_pos,effect_allele,other_allele,beta_qtl,se_qtl,p_qtl,gwas_beta_h,gwas_se,gwas_p,relation)],
      file.path(cout,"harmonized_locus.tsv"),sep="\t"
    )

    rec<-data.table(
      tissue,circ_id=circ,
      circ_chr=cm$circ_chr_key,circ_start=cm$circ_start_key,circ_end=cm$circ_end_key,
      status=ans$status,n_locus=n_locus,n_shared=nrow(x),shared_fraction=nrow(x)/n_locus,
      N_cQTL=N_CQTL,N_GWAS=GWAS_N,
      n_candidate_SNPs=cm$n_candidate_SNPs,candidate_SNPs=cm$candidate_SNPs,
      best_candidate_SNP=cm$best_candidate_SNP,best_ALS_p=cm$best_ALS_p,
      best_ALS_tissue_FDR=cm$best_ALS_tissue_FDR,best_ALS_global_FDR=cm$best_ALS_global_FDR,
      cqtl_LD_mismatch_s=s_cqtl,gwas_LD_mismatch_s=s_gwas,
      n_cQTL_CS=ans$ncs1,n_GWAS_CS=ans$ncs2,
      top_cQTL_SNP=ans$cqtl_top,top_cQTL_PIP=ans$cqtl_pip,
      top_GWAS_SNP=ans$gwas_top,top_GWAS_PIP=ans$gwas_pip,
      max_PP_H4=ans$h4,H4_reference_pass=!is.na(ans$h4)&ans$h4>=H4REF
    )

    master[[length(master)+1]]<-rec

    cat(
      ans$status,", shared=",nrow(x),"/",n_locus,
      " (",round(100*nrow(x)/n_locus,1),"%), candidate_SNPs=",cm$n_candidate_SNPs,
      ", cQTL_CS=",ans$ncs1,", GWAS_CS=",ans$ncs2,", H4=",signif(ans$h4,3),"\n",sep=""
    )

    rm(G,R,x,q); gc()
  }

  rm(full,cloc,cand,circ_meta); gc()
}

# ==============================================================
# MAIN SUMMARY
# ==============================================================
res<-if(length(master)) rbindlist(master,fill=TRUE) else data.table()
if(!"max_PP_H4"%in%names(res)) res[,max_PP_H4:=NA_real_]

outfile<-file.path(BASE,"MR","cQTL_SuSiE_coloc_summary.tsv")
fwrite(res,outfile,sep="\t")

# ==============================================================
# STATUS COUNTS
# ==============================================================
status_counts<-if(nrow(res)){
  z<-copy(res)
  z[,status_group:=ifelse(grepl("^ERROR:",status),"ERROR",status)]
  z[,.(n_loci=.N),by=.(status=status_group)][order(-n_loci)]
} else data.table(status=character(),n_loci=integer())

if(nrow(status_counts)) status_counts[,percent:=round(100*n_loci/sum(n_loci),2)] else status_counts[,percent:=numeric()]

statusfile<-file.path(BASE,"MR","cQTL_SuSiE_coloc_status_counts.tsv")
fwrite(status_counts,statusfile,sep="\t")

# ==============================================================
# PP.H4 VALUES + HISTOGRAM
# ==============================================================
h4<-if(nrow(res)){
  res[is.finite(max_PP_H4),.(tissue,circ_id,max_PP_H4,n_candidate_SNPs,best_candidate_SNP,best_ALS_p,best_ALS_global_FDR)]
} else data.table()

h4file<-file.path(BASE,"MR","cQTL_PPH4_values.tsv")
plotfile<-file.path(BASE,"MR","cQTL_PPH4_histogram.png")
fwrite(h4,h4file,sep="\t")

if(nrow(h4)){
  p<-ggplot(h4,aes(x=max_PP_H4))+
    geom_histogram(binwidth=.05,boundary=0,closed="left")+
    geom_vline(xintercept=H4REF,linetype="dashed",linewidth=.8)+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.1))+
    labs(
      title=paste0("cQTL SuSiE colocalization PP.H4 distribution (n=",nrow(h4),")"),
      subtitle=paste0("Dashed line = ",H4REF," reference threshold"),
      x="Maximum PP.H4 per circRNA-tissue locus",
      y="Number of loci"
    )+
    theme_bw(base_size=12)

  ggsave(filename=plotfile,plot=p,width=8,height=5,dpi=300)

  cat(
    "\nPP.H4 DISTRIBUTION\n",
    "n = ",nrow(h4),
    "\nmedian = ",median(h4$max_PP_H4),
    "\nIQR = ",quantile(h4$max_PP_H4,.25)," to ",quantile(h4$max_PP_H4,.75),
    "\n>=0.50 = ",sum(h4$max_PP_H4>=.50),
    "\n>=0.80 = ",sum(h4$max_PP_H4>=.80),
    "\nValues: ",h4file,
    "\nHistogram: ",plotfile,"\n",
    sep=""
  )
} else {
  cat("\nNo finite PP.H4 values; histogram not generated.\n")
}

cat("\n",strrep("=",70),"\nSTATUS COUNTS\n",strrep("=",70),"\n",sep="")
print(status_counts)

cat("\n",strrep("=",70),"\ncQTL SuSiE / COLOC SUMMARY\n",strrep("=",70),"\n",sep="")
print(res)

cat(
  "\nMain summary: ",outfile,
  "\nStatus counts: ",statusfile,
  "\nPP.H4 values: ",h4file,"\n",
  sep=""
)

if(nrow(h4)) cat("Histogram: ",plotfile,"\n",sep="")
RS

echo
echo "cQTL SuSiE + colocalization complete."
echo
echo "=== STATUS COUNTS ==="
column -t -s $'\t' "$BASE/MR/cQTL_SuSiE_coloc_status_counts.tsv"

echo
echo "=== OUTPUT FILES ==="
echo "$BASE/MR/cQTL_SuSiE_coloc_summary.tsv"
echo "$BASE/MR/cQTL_SuSiE_coloc_status_counts.tsv"
echo "$BASE/MR/cQTL_PPH4_values.tsv"
echo "$BASE/MR/cQTL_PPH4_histogram.png"
