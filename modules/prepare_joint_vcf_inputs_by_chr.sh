#!/bin/bash
# AnswerALS parallel, chromosome-by-chromosome replacement for prepare_joint_vcf_inputs.sh (takes too long for gVCFs)
#
# Run this from a login node with:
#   bash prepare_joint_vcf_inputs_by_chr.sh
#
# It submits one preparation job, a 25-task chromosome array, and a finalizer.

set -euo pipefail
shopt -s nullglob

EXPERIMENT="answer_ALS"
BASE="/home/zw529/donglab/data/${EXPERIMENT}"
QTL_DIR="$BASE/QTL"
RNA_BASE="$BASE/iPSC_Motor_Neuron/RNAseq/Processed"
WGS_BASE="$BASE/Blood/WGS/Processed"
REFERENCE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
EXPECTED_GVCFS=815
THREADS=8
MAX_CONCURRENT_CHR_JOBS=4
FINAL_VCF="$QTL_DIR/joint_genotyped_GQ.vcf.gz"
CHROM_DIR="$QTL_DIR/chromosome_joint_vcfs"

CHROMS=(
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
)

load_modules() {
    module load BCFtools/1.21-GCC-13.3.0
    module load HTSlib/1.21-GCC-13.3.0
}

prepare_inputs() {
    local vcf_list="$RUN_DIR/vcf_merge_list.txt"
    local sample_map="$RUN_DIR/sample_map.txt"
    local chrom_list="$RUN_DIR/chromosomes.txt"
    local gvcf_count=0
    local contig_signature=""
    local donor source_gvcf vcf_id current_contigs header sample_count
    declare -A seen_vcf_ids=()

    [[ -s "$REFERENCE" ]] || {
        echo "[ERROR] Reference FASTA not found: $REFERENCE" >&2
        exit 1
    }
    [[ -s "${REFERENCE}.fai" ]] || {
        echo "[ERROR] Reference FASTA index not found: ${REFERENCE}.fai" >&2
        exit 1
    }

    for chr in "${CHROMS[@]}"; do
        if ! cut -f1 "${REFERENCE}.fai" | grep -Fxq "$chr"; then
            echo "[ERROR] Reference index is missing $chr: ${REFERENCE}.fai" >&2
            exit 1
        fi
    done

    printf '%s\n' "${CHROMS[@]}" > "$chrom_list"
    : > "$vcf_list"
    : > "$sample_map"

    echo "Preparing AnswerALS raw gVCFs..."
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
        vcf_id=$(bcftools query -l "$source_gvcf")
        sample_count=$(printf '%s\n' "$vcf_id" | awk 'NF {n++} END {print n+0}')

        if ((sample_count != 1)); then
            echo "[ERROR] Expected one sample in $source_gvcf; found $sample_count." >&2
            exit 1
        fi

        if [[ -n "${seen_vcf_ids[$vcf_id]+x}" ]]; then
            echo "[ERROR] Duplicate VCF sample ID '$vcf_id' for $donor and ${seen_vcf_ids[$vcf_id]}." >&2
            exit 1
        fi
        seen_vcf_ids["$vcf_id"]="$donor"

        current_contigs=$(bcftools view -h "$source_gvcf" |
            sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p')

        if [[ -z "$current_contigs" ]]; then
            echo "[ERROR] No contig declarations found in $source_gvcf" >&2
            exit 1
        fi

        if ((gvcf_count == 0)); then
            contig_signature="$current_contigs"
            header=$(bcftools view -h "$source_gvcf")

            if [[ "$header" != *"##ALT=<ID=NON_REF"* ||
                  "$header" != *"##INFO=<ID=END"* ]]; then
                echo "[ERROR] Missing expected gVCF NON_REF or INFO/END header in $source_gvcf" >&2
                exit 1
            fi

            for chr in "${CHROMS[@]}"; do
                if ! grep -Fxq "$chr" <<< "$current_contigs"; then
                    echo "[ERROR] $source_gvcf does not declare $chr." >&2
                    exit 1
                fi
            done
        elif [[ "$current_contigs" != "$contig_signature" ]]; then
            echo "[ERROR] Contig order differs in $source_gvcf." >&2
            exit 1
        fi

        printf '%s\n' "$source_gvcf" >> "$vcf_list"
        printf '%s\t%s\n' "$vcf_id" "$donor" >> "$sample_map"

        if [[ ! -s "${source_gvcf}.csi" && ! -s "${source_gvcf}.tbi" ]]; then
            bcftools index --csi --threads "$THREADS" "$source_gvcf"
        fi

        gvcf_count=$((gvcf_count + 1))
    done

    if ((gvcf_count != EXPECTED_GVCFS)); then
        echo "[ERROR] Expected $EXPECTED_GVCFS raw gVCFs; found $gvcf_count. Refusing incomplete cohort." >&2
        exit 1
    fi

    echo "Prepared $gvcf_count raw gVCFs."
}

merge_chromosome() {
    local task_id="${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
    local chrom="${CHROMS[$task_id]}"
    local vcf_list="$RUN_DIR/vcf_merge_list.txt"
    local sample_map="$RUN_DIR/sample_map.txt"
    local merged_tmp="$RUN_DIR/.${EXPERIMENT}_${chrom}.${SLURM_JOB_ID}.merge.tmp.vcf.gz"
    local shard_tmp="$RUN_DIR/.${EXPERIMENT}_${chrom}.${SLURM_JOB_ID}.reheader.tmp.vcf.gz"
    local shard="$RUN_DIR/${EXPERIMENT}_${chrom}.vcf.gz"

    [[ -s "$vcf_list" && -s "$sample_map" ]] || {
        echo "[ERROR] Prepared input lists are missing in $RUN_DIR." >&2
        exit 1
    }

    trap 'rm -f "$merged_tmp" "$shard_tmp" "${shard_tmp}.tbi" "${shard_tmp}.csi"' EXIT

    echo "Merging $chrom from $EXPECTED_GVCFS raw gVCFs..."
    bcftools merge \
        --gvcf "$REFERENCE" \
        --regions "$chrom" \
        --merge 'both,*' \
        --file-list "$vcf_list" \
        --output-type u |
    bcftools view \
        --types snps,indels,mnps \
        --output-type u |
    bcftools annotate \
        -x INFO,^FORMAT/GT,FORMAT/GQ \
        --threads "$THREADS" \
        --output-type z \
        --output "$merged_tmp"

    echo "Reheadering $chrom to canonical contigs and AnswerALS donor IDs..."
    bcftools reheader \
        --fai "${REFERENCE}.fai" \
        --samples "$sample_map" \
        --threads "$THREADS" \
        --output "$shard_tmp" \
        "$merged_tmp"

    bcftools index -f -t "$shard_tmp"
    mv "$shard_tmp" "$shard"
    mv "${shard_tmp}.tbi" "${shard}.tbi"
    rm -f "$merged_tmp"
    trap - EXIT

    echo "Finished $chrom: $shard"
}

finalize_joint_vcf() {
    local shard_list="$RUN_DIR/chromosome_vcf_list.txt"
    local joint_tmp="$QTL_DIR/joint_genotyped_GQ.${RUN_ID}.tmp.vcf.gz"
    local sample_count variant_count

    : > "$shard_list"
    for chr in "${CHROMS[@]}"; do
        local shard="$RUN_DIR/${EXPERIMENT}_${chr}.vcf.gz"
        [[ -s "$shard" ]] || {
            echo "[ERROR] Missing chromosome shard: $shard" >&2
            exit 1
        }
        [[ -s "${shard}.tbi" || -s "${shard}.csi" ]] || {
            echo "[ERROR] Missing index for chromosome shard: $shard" >&2
            exit 1
        }
        printf '%s\n' "$shard" >> "$shard_list"
    done

    trap 'rm -f "$joint_tmp" "${joint_tmp}.tbi" "${joint_tmp}.csi"' EXIT

    echo "Concatenating chromosome VCFs..."
    bcftools concat \
        --file-list "$shard_list" \
        --threads "$THREADS" \
        --output-type z \
        --output "$joint_tmp"

    sample_count=$(bcftools query -l "$joint_tmp" | awk 'END {print NR}')
    if ((sample_count != EXPECTED_GVCFS)); then
        echo "[ERROR] Concatenated VCF has $sample_count samples; expected $EXPECTED_GVCFS." >&2
        exit 1
    fi

    bcftools index -f -t "$joint_tmp"
    variant_count=$(bcftools index --stats "$joint_tmp" |
        awk '{n += $3} END {print n+0}')

    if ((variant_count == 0)); then
        echo "[ERROR] Concatenated VCF contains no variant records." >&2
        exit 1
    fi

    echo "Installing final joint VCF..."
    mv -f "$joint_tmp" "$FINAL_VCF"
    mv -f "${joint_tmp}.tbi" "${FINAL_VCF}.tbi"
    rm -f "${FINAL_VCF}.csi"
    trap - EXIT

    echo "AnswerALS samples in final VCF: $sample_count"
    echo "AnswerALS variant records: $variant_count"
    echo "Chromosome VCFs retained in: $RUN_DIR"
    echo "Final joint VCF: $FINAL_VCF"
}

MODE="${1:-submit}"

if [[ "$MODE" == "submit" ]]; then

    mkdir -p "$CHROM_DIR"
    RUN_ID="answerALS_$(date +%Y%m%d_%H%M%S)_$$"
    RUN_DIR="$CHROM_DIR/$RUN_ID"
    mkdir -p "$RUN_DIR/logs"
    SCRIPT_PATH=$(readlink -f "$0")

    PREP_JOB=$(sbatch --parsable \
        --job-name=answerALS_gvcf_prep \
        --partition=day \
        --time=23:00:00 \
        --cpus-per-task=4 \
        --mem=16G \
        --output="$RUN_DIR/logs/prepare_%j.out" \
        --error="$RUN_DIR/logs/prepare_%j.err" \
        --export=ALL,BYCHR_RUN_ID="$RUN_ID" \
        "$SCRIPT_PATH" _prepare)
    PREP_JOB="${PREP_JOB%%;*}"

    ARRAY_JOB=$(sbatch --parsable \
        --job-name=answerALS_gvcf_chr \
        --partition=week \
        --time=5-23:00:00 \
        --cpus-per-task=8 \
        --mem=99G \
        --array="0-$((${#CHROMS[@]} - 1))%${MAX_CONCURRENT_CHR_JOBS}" \
        --dependency="afterok:${PREP_JOB}" \
        --output="$RUN_DIR/logs/chr_%A_%a.out" \
        --error="$RUN_DIR/logs/chr_%A_%a.err" \
        --export=ALL,BYCHR_RUN_ID="$RUN_ID" \
        "$SCRIPT_PATH" _worker)
    ARRAY_JOB="${ARRAY_JOB%%;*}"

    FINAL_JOB=$(sbatch --parsable \
        --job-name=answerALS_gvcf_finalize \
        --partition=week \
        --time=4-00:00:00 \
        --cpus-per-task=8 \
        --mem=16G \
        --dependency="afterok:${ARRAY_JOB}" \
        --output="$RUN_DIR/logs/finalize_%j.out" \
        --error="$RUN_DIR/logs/finalize_%j.err" \
        --export=ALL,BYCHR_RUN_ID="$RUN_ID" \
        "$SCRIPT_PATH" _finalize)
    FINAL_JOB="${FINAL_JOB%%;*}"

    echo "Submitted AnswerALS gVCF workflow."
    echo "Preparation job: $PREP_JOB"
    echo "Chromosome array: $ARRAY_JOB (${#CHROMS[@]} tasks; at most ${MAX_CONCURRENT_CHR_JOBS} running at once)"
    echo "Finalization job: $FINAL_JOB (runs only if every chromosome task succeeds)"
    echo "Chromosome VCFs: $RUN_DIR"
    echo "Final joint VCF: $FINAL_VCF"
    exit 0
fi

RUN_ID="${BYCHR_RUN_ID:?Missing BYCHR_RUN_ID from Slurm submission}"
RUN_DIR="$CHROM_DIR/$RUN_ID"

load_modules

case "$MODE" in
    _prepare)
        mkdir -p "$RUN_DIR"
        prepare_inputs
        ;;
    _worker)
        merge_chromosome
        ;;
    _finalize)
        finalize_joint_vcf
        ;;
    *)
        echo "[ERROR] Unknown mode: $MODE" >&2
        exit 1
        ;;
esac
