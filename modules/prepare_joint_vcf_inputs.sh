#!/bin/bash
#SBATCH --job-name=prep_joint_vcf
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/prep_joint_vcf.out   ### change log paths based on experiment
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/prep_joint_vcf.err   ### change log paths based on experiment
#SBATCH --time=5-24:00:00
#SBATCH --partition=week
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

set -euo pipefail
shopt -s nullglob

module load Python/3.12.3-GCCcore-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0

# Set to target_ALS or answer_ALS before running.
EXPERIMENT="answer_ALS"

BASE="/home/zw529/donglab/data/$EXPERIMENT"
QTL_DIR="$BASE/QTL"
mkdir -p "$QTL_DIR"

VCF_LIST="$QTL_DIR/vcf_merge_list.txt"
SAMPLE_MAP="$QTL_DIR/sample_map.txt"
FINAL_VCF="$QTL_DIR/joint_genotyped_GQ.vcf.gz"

# --- STEP 1: SELECT VCFs AND PREPARE INPUTS ---

if [[ "$EXPERIMENT" == "target_ALS" ]]; then

    # Keep the existing TargetALS VCF ranking and sidecar behavior.
    echo "Ranking TargetALS VCFs and creating standardized sidecars..."

    find "$BASE" -name "*.vcf.gz" -printf '%h\n' | sort -u | while read -r dir; do
        BEST_VCF=$(find "$dir" -maxdepth 1 -name "*.vcf.gz" \
            ! -path "*/QTL/*" \
            ! -name "*annotated*" \
            ! -name "*cnv*" \
            ! -name "*sv*" \
            ! -name "*repeats*" \
            ! -name "*manta*" \
            ! -name "*canvas*" \
            ! -name "*.g.vcf.gz" \
            ! -name "genotypes_quality.vcf.gz" | \
            xargs -r ls -S | \
            awk '
                /recalibrated/ { print 1, $0; next }
                /hard-filtered/ { print 2, $0; next }
                { print 3, $0 }
            ' | sort -n | head -n 1 | cut -d' ' -f2-)

        if [[ -n "$BEST_VCF" ]]; then
            SIDECAR="$dir/genotypes_quality.vcf.gz"
            if [[ ! -f "$SIDECAR" ]]; then
                echo "Processing winner in $dir: $(basename "$BEST_VCF")"
                bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ \
                    "$BEST_VCF" -O z -o "$SIDECAR"
                tabix -f -p vcf "$SIDECAR"
            fi
        fi
    done

    # --- STEP 2: TARGETALS METADATA RECONCILIATION ---
    python3 - <<'EOF'
import pandas as pd
from pathlib import Path
import subprocess

BASE = Path("/home/zw529/donglab/data/target_ALS")
QTL_DIR = BASE / "QTL"

wgs_meta = pd.read_csv(BASE / "targetALS_wgs_metadata.csv", encoding="utf-8-sig")
rnaseq_meta = pd.read_csv(BASE / "targetALS_rnaseq_metadata.csv")
wgs_meta.columns = wgs_meta.columns.str.strip().str.lower()
rnaseq_meta.columns = rnaseq_meta.columns.str.strip().str.lower()

shared = set(wgs_meta["externalsubjectid"].dropna()) & \
         set(rnaseq_meta["externalsubjectid"].dropna())
wgs = wgs_meta[wgs_meta["externalsubjectid"].isin(shared)].copy()
sidecars = list(BASE.glob("**/genotypes_quality.vcf.gz"))

def find_sidecar(sample_id):
    sid = str(sample_id).lower()
    sid_alt = sid.replace("-", "_")
    for path in sidecars:
        p = str(path).lower()
        if path.name == "genotypes_quality.vcf.gz" and (sid in p or sid_alt in p):
            return str(path)
    return None

wgs["vcf_path"] = wgs["externalsampleid"].apply(find_sidecar)
found = wgs[wgs["vcf_path"].notna()].drop_duplicates("externalsubjectid")

print(f"Final subject count for merge: {len(found)}")
found["vcf_path"].to_csv(QTL_DIR / "vcf_merge_list.txt", index=False, header=False)

sample_map = []
for _, row in found.iterrows():
    result = subprocess.run(
        ["bcftools", "query", "-l", row["vcf_path"]],
        capture_output=True, text=True, check=True
    )
    vcf_id = result.stdout.strip()
    if vcf_id:
        sample_map.append({"vcf_id": vcf_id, "subject_id": row["externalsubjectid"]})

pd.DataFrame(sample_map).to_csv(
    QTL_DIR / "sample_map.txt", sep="\t", index=False, header=False
)
EOF

    # --- STEP 3: TARGETALS MERGE & REHEADER ---
    echo "Merging TargetALS sidecars (preserving GT and GQ)..."
    bcftools merge -f PASS -l "$VCF_LIST" \
        --threads 8 -O z -o "$FINAL_VCF"

    echo "Reheadering sequencing IDs to subject IDs..."
    bcftools reheader -s "$SAMPLE_MAP" \
        -o "${FINAL_VCF}.tmp" "$FINAL_VCF"
    mv "${FINAL_VCF}.tmp" "$FINAL_VCF"
    tabix -f -p vcf "$FINAL_VCF"

elif [[ "$EXPERIMENT" == "answer_ALS" ]]; then

    echo "Preparing AnswerALS raw gVCFs for joint genotyping..."

    module load GATK/4.6.1.0-GCCcore-13.3.0-Java-17

    RNA_BASE="$BASE/iPSC_Motor_Neuron/RNAseq/Processed"
    WGS_BASE="$BASE/Blood/WGS/Processed"
    REFERENCE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
    EXPECTED_GVCFS=815

    CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

    [[ -s "$REFERENCE" ]] || {
        echo "[ERROR] Reference FASTA not found: $REFERENCE" >&2
        exit 1
    }
    [[ -s "${REFERENCE}.fai" ]] || {
        echo "[ERROR] Reference FASTA index not found: ${REFERENCE}.fai" >&2
        exit 1
    }
    [[ -s "${REFERENCE%.fa}.dict" ]] || {
        echo "[ERROR] GATK sequence dictionary not found: ${REFERENCE%.fa}.dict" >&2
        exit 1
    }

    WORK_DIR="$QTL_DIR/gatk_joint_work_${SLURM_JOB_ID:-manual}"
    CLEAN_GVCF_DIR="$WORK_DIR/gvcfs"
    DB_DIR="$WORK_DIR/genomicsdb"
    RAW_JOINT_VCF="$WORK_DIR/joint_genotyped.raw.vcf.gz"
    GATK_SAMPLE_MAP="$WORK_DIR/gatk_sample_name_map.tsv"

    mkdir -p "$CLEAN_GVCF_DIR" "$WORK_DIR/tmp"
    : > "$GATK_SAMPLE_MAP"
    : > "$SAMPLE_MAP"

    declare -A seen_vcf_ids
    gvcf_count=0

    for rna_dir in "$RNA_BASE"/*/; do
        [[ -d "$rna_dir" ]] || continue
        donor=$(basename "$rna_dir")
        candidates=("$WGS_BASE/$donor"/*.haplotypeCalls.er.raw.vcf.gz)

        if ((${#candidates[@]} == 0)); then
            echo "[INFO] $donor: no raw gVCF in its WGS folder; skipping."
            continue
        elif ((${#candidates[@]} != 1)); then
            echo "[ERROR] $donor: expected one raw gVCF; found ${#candidates[@]}." >&2
            exit 1
        fi

        source_gvcf="${candidates[0]}"
        vcf_samples=$(bcftools query -l "$source_gvcf")
        sample_count=$(printf '%s\n' "$vcf_samples" | awk 'NF {n++} END {print n+0}')

        if ((sample_count != 1)); then
            echo "[ERROR] $donor: expected one VCF sample, found $sample_count in $source_gvcf" >&2
            exit 1
        fi

        vcf_id="$vcf_samples"
        if [[ -n "${seen_vcf_ids[$vcf_id]+x}" ]]; then
            echo "[ERROR] Duplicate VCF sample ID '$vcf_id' found for $donor and ${seen_vcf_ids[$vcf_id]}." >&2
            exit 1
        fi
        seen_vcf_ids["$vcf_id"]="$donor"

        clean_gvcf="$CLEAN_GVCF_DIR/${donor}.g.vcf.gz"
        subset_gvcf="$CLEAN_GVCF_DIR/${donor}.subset.tmp.vcf.gz"

        echo "Preparing $donor from $(basename "$source_gvcf")"

        # Keep only the standard chromosomes requested for this QTL workflow.
        # -t streams the source gVCF, so a source index is not required.
        bcftools view -t "$CHROMS" -Oz \
            -o "$subset_gvcf" "$source_gvcf"

        # Match the contig declarations to the selected reference FASTA.
        # This removes header contigs absent from the FASTA index.
        bcftools reheader -f "${REFERENCE}.fai" \
            -o "$clean_gvcf" "$subset_gvcf"
        rm "$subset_gvcf"
        bcftools index -f -t "$clean_gvcf"

        printf '%s\t%s\n' "$vcf_id" "$clean_gvcf" >> "$GATK_SAMPLE_MAP"
        printf '%s\t%s\n' "$vcf_id" "$donor" >> "$SAMPLE_MAP"
        ((gvcf_count+=1))
    done

    if ((gvcf_count != EXPECTED_GVCFS)); then
        echo "[ERROR] Expected $EXPECTED_GVCFS AnswerALS raw gVCFs, found $gvcf_count. Refusing to run an incomplete cohort." >&2
        exit 1
    fi

    echo "Prepared $gvcf_count raw gVCFs."
    echo "Running GATK GenomicsDBImport..."

    interval_args=()
    IFS=',' read -r -a chrom_array <<< "$CHROMS"
    for chr in "${chrom_array[@]}"; do
        interval_args+=(-L "$chr")
    done

    gatk --java-options "-Xms32g -Xmx32g" GenomicsDBImport \
        -R "$REFERENCE" \
        --sample-name-map "$GATK_SAMPLE_MAP" \
        --genomicsdb-workspace-path "$DB_DIR" \
        --batch-size 50 \
        --reader-threads 4 \
        --tmp-dir "$WORK_DIR/tmp" \
        "${interval_args[@]}"

    echo "Running GATK GenotypeGVCFs..."
    gatk --java-options "-Xms32g -Xmx32g" GenotypeGVCFs \
        -R "$REFERENCE" \
        -V "gendb://${DB_DIR}" \
        -O "$RAW_JOINT_VCF" \
        "${interval_args[@]}"

    echo "Keeping GT and GQ fields and renaming samples to AnswerALS donor IDs..."
    bcftools annotate -x INFO,^FORMAT/GT,FORMAT/GQ \
        "$RAW_JOINT_VCF" -Oz -o "${FINAL_VCF}.tmp"

    bcftools reheader -s "$SAMPLE_MAP" \
        -o "$FINAL_VCF" "${FINAL_VCF}.tmp"
    rm "${FINAL_VCF}.tmp"
    tabix -f -p vcf "$FINAL_VCF"

    final_sample_count=$(bcftools query -l "$FINAL_VCF" | wc -l)
    if ((final_sample_count != EXPECTED_GVCFS)); then
        echo "[ERROR] Final VCF has $final_sample_count samples; expected $EXPECTED_GVCFS." >&2
        exit 1
    fi

else
    echo "[ERROR] EXPERIMENT must be target_ALS or answer_ALS." >&2
    exit 1
fi

echo "Processing complete. Joint VCF ready at $FINAL_VCF"
echo "Sample count: $(bcftools query -l "$FINAL_VCF" | wc -l)"
