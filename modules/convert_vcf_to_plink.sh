#!/bin/bash
#SBATCH --job-name=vcf_to_plink_imputed
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/vcf_to_plink.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/vcf_to_plink.err

set -euo pipefail

#----------------------------------------
# PATHS & PREFIXES
#----------------------------------------
# VCF_IN=/home/zw529/donglab/data/target_ALS/QTL/target_ALS_imputed_filtered_joint.vcf.gz
# OUTDIR=/home/zw529/donglab/data/target_ALS/QTL/plink
VCF_IN=/home/zw529/donglab/data/answer_ALS/QTL/joint_genotyped_GQ.vcf.gz   ### imputation not performed on answerALS .vcf; too dense
OUTDIR=/home/zw529/donglab/data/answer_ALS/QTL/plink
mkdir -p ${OUTDIR}

# AnswerALS sex comes from the participant metadata
SEX_PSAM=${OUTDIR}/sexcheck/answerALS_sex.psam
ANSWERALS_MODE=false
if [[ "$VCF_IN" == *"/answer_ALS/"* ]]; then
    ANSWERALS_MODE=true
fi

#----------------------------------------
# OPTIONAL ANSWERALS SEX-CHECK FUNCTION
# Uncomment this function and its call to recreate the AnswerALS metadata PSAM and chrX/chrY sex-check.
#----------------------------------------
# answerals_sexcheck() {
#     local W="${OUTDIR}/sexcheck"
#     local V="${VCF_IN}"
#     local M="/home/zw529/donglab/data/answer_ALS/aals_dataportal_datatable.csv"
#
#     mkdir -p "$W"
#     module load BCFtools/1.21-GCC-13.3.0
#
#     python3 -c 'import csv,subprocess,pathlib,sys; m=pathlib.Path(sys.argv[1]); v=sys.argv[2]; p=pathlib.Path(sys.argv[3]); r=csv.DictReader(m.open(encoding="utf-8-sig",newline="")); idc=r.fieldnames[0]; sc=next(c for c in r.fieldnames if c.strip().lower()=="sex"); d={x[idc].strip():x[sc].strip().lower() for x in r if x.get(idc)}; ids=subprocess.check_output(["bcftools","query","-l",v],text=True).splitlines(); assert all(i in d for i in ids),"VCF sample missing from metadata"; p.write_text("#IID SEX\n"+"".join("{} {}\n".format(i,{"male":1,"female":2}.get(d[i],0)) for i in ids))' \
#         "/home/zw529/donglab/data/answer_ALS/aals_dataportal_datatable.csv" \
#         "$V" "$W/answerALS_sex.psam"
#
#     bcftools view -r chrX,chrY -Oz -o "$W/sex_chroms.vcf.gz" "$V"
#     bcftools index -t "$W/sex_chroms.vcf.gz"
#     module load PLINK2/avx2_20250707
#     plink2 --vcf "$W/sex_chroms.vcf.gz" \
#         --psam "$W/answerALS_sex.psam" \
#         --split-par hg38 \
#         --check-sex max-female-xf=0.2 min-male-xf=0.8 \
#         --out "$W/answerALS_sexcheck"
# }
# answerals_sexcheck

# Temp files for Step 0 (if ever uncommented)
VCF_FILTERED=${OUTDIR}/joint_GQ0_temp.vcf.gz
VCF_BIALLELIC=${OUTDIR}/joint_biallelic_temp.vcf.gz

# --- OLD AUTOSOME-ONLY PREFIXES (HASHED OUT) ---
# RAW_PREFIX=${OUTDIR}/joint_autosomes_raw
# QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
# FILTERED_PREFIX=${OUTDIR}/joint_autosomes_filtered
# BED_PREFIX=${OUTDIR}/joint_autosomes_filtered_bed
# MATRIX_PREFIX=${OUTDIR}/joint_autosomes_matrixEQTL

# --- NEW ALL-CHROMOSOME PREFIXES ---
RAW_PREFIX=${OUTDIR}/joint_all_chrs_raw
QC_SAMPLE_PREFIX=${OUTDIR}/joint_samples_qc
FILTERED_PREFIX=${OUTDIR}/joint_all_chrs_filtered
BED_PREFIX=${OUTDIR}/joint_all_chrs_filtered_bed
MATRIX_PREFIX=${OUTDIR}/joint_all_chrs_matrixEQTL

#----------------------------------------
# STEP 0: PRE-PROCESSING (OBSOLETE FOR IMPUTED VCFs)
#----------------------------------------
# Only uncomment if switching back to "Normal" raw VCFs.
#
# echo "Running Step 0a (GQ Filter)..."
# module --force purge
# module load StdEnv
# module load BCFtools/1.21
# bcftools filter -e 'FORMAT/GQ < 0' -S . ${VCF_IN} -O z -o ${VCF_FILTERED}
# bcftools index -t ${VCF_FILTERED}
#
# echo "Running Step 0b (Multiallelic to Biallelic)..."
# bcftools norm -m -any ${VCF_FILTERED} -O z -o ${VCF_BIALLELIC}
# bcftools index -t ${VCF_BIALLELIC}
# rm ${VCF_FILTERED} ${VCF_FILTERED}.tbi
# VCF_IN=${VCF_BIALLELIC}

#----------------------------------------
# STEP 1: VCF → PLINK2 IMPORT
#----------------------------------------
module --force purge
module load StdEnv
module load PLINK2/avx2_20250707

# Keep this option documented. TargetALS still uses it below;
# AnswerALS uses the metadata sex values supplied in SEX_PSAM.
#     --impute-sex max-female-xf=0.2 min-male-xf=0.8 \

#================================================================================
# VERSION 1: WITH SEX CHROMOSOMES (ACTIVE)
#================================================================================
echo "Importing VCF with Sex Chromosomes..."

if "$ANSWERALS_MODE"; then
    [[ -s "$SEX_PSAM" ]] || {
        echo "[ERROR] AnswerALS sex PSAM is missing or empty: $SEX_PSAM" >&2
        echo "[ERROR] Create it with the optional answerals_sexcheck function above." >&2
        exit 1
    }

    plink2 --vcf ${VCF_IN} \
        --psam ${SEX_PSAM} \
        --chr chr1-22,chrX,chrY,chrM \
        --split-par hg38 \
        --set-all-var-ids @:#:\$r:\$a \
        --new-id-max-allele-len 800 \
        --vcf-half-call m \
        --import-max-alleles 255 \
        --make-pgen \
        --out ${RAW_PREFIX} \
        --threads ${SLURM_CPUS_PER_TASK}
else
    plink2 --vcf ${VCF_IN} \
        --chr chr1-22,chrX,chrY,chrM \
        --split-par hg38 \
        --set-all-var-ids @:#:\$r:\$a \
        --new-id-max-allele-len 800 \
        --vcf-half-call m \
        --impute-sex max-female-xf=0.2 min-male-xf=0.8 \
        --make-pgen \
        --out ${RAW_PREFIX} \
        --threads ${SLURM_CPUS_PER_TASK}
fi

echo "Running sex-check validation..."
plink2 --pfile ${RAW_PREFIX} --check-sex max-female-xf=0.2 min-male-xf=0.8 --out ${QC_SAMPLE_PREFIX}

if "$ANSWERALS_MODE"; then
    # Keep sex-check results for review; ambiguous results are not automatic exclusions.
    : > ${OUTDIR}/fail_sex.txt
else
    grep "PROBLEM" ${QC_SAMPLE_PREFIX}.sexcheck | awk '{print $1, $2}' > ${OUTDIR}/fail_sex.txt || touch ${OUTDIR}/fail_sex.txt
fi


#================================================================================
# VERSION 2: WITHOUT SEX CHROMOSOMES (HASHED OUT)
#================================================================================
# echo "Importing Autosome-only VCF..."
# plink2 --vcf ${VCF_IN} \
#        --chr chr1-22 \
#        --set-all-var-ids @:#:\$r:\$a \
#        --new-id-max-allele-len 800 \
#        --vcf-half-call m \
#        --make-pgen \
#        --out ${RAW_PREFIX} \
#        --threads ${SLURM_CPUS_PER_TASK}
#
# touch ${OUTDIR}/fail_sex.txt

#----------------------------------------
# SAMPLE-LEVEL QC (Shared Logic)
#----------------------------------------
echo "Performing sample-level checks..."
plink2 --pfile ${RAW_PREFIX} --het --out ${QC_SAMPLE_PREFIX}
plink2 --pfile ${RAW_PREFIX} --king-cutoff 0.45 --out ${QC_SAMPLE_PREFIX}_relatedness

awk 'NR>1 && ($6 > 0.2 || $6 < -0.2) {print $1, $2}' ${QC_SAMPLE_PREFIX}.het > ${OUTDIR}/fail_het.txt || touch ${OUTDIR}/fail_het.txt

if [ -f ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id ]; then
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt ${QC_SAMPLE_PREFIX}_relatedness.king.cutoff.out.id | sort | uniq > ${OUTDIR}/all_fail_samples.txt
else
    cat ${OUTDIR}/fail_sex.txt ${OUTDIR}/fail_het.txt | sort | uniq > ${OUTDIR}/all_fail_samples.txt
fi

#----------------------------------------
# FINAL QC FILTERING
#----------------------------------------
echo "Applying variant missingness filter (--geno 0.05)..."
plink2 --pfile ${RAW_PREFIX} \
       --geno 0.05 \
       --make-pgen \
       --out ${OUTDIR}/temp_variant_filtered

echo "Applying final hard sample filters, MAF, and HWE thresholds..."
plink2 --pfile ${OUTDIR}/temp_variant_filtered \
       --remove ${OUTDIR}/all_fail_samples.txt \
       --mind 0.05 \
       --maf 0.05 \
       --hwe 1e-6 \
       --max-alleles 2 \
       --make-pgen \
       --out ${FILTERED_PREFIX}

#----------------------------------------
# PCA & OUTPUT
#----------------------------------------
echo "Running PCA and generating final formats..."
plink2 --pfile ${FILTERED_PREFIX} --pca 10 --out ${OUTDIR}/joint_pca
plink2 --pfile ${FILTERED_PREFIX} --make-bed --out ${BED_PREFIX}

plink2 --pfile ${FILTERED_PREFIX} \
       --recode A \
       --hard-call-threshold 0.49 \
       --nonfounders \
       --out ${MATRIX_PREFIX}

echo "Pipeline Finished."
