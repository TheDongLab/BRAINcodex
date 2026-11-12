#!/usr/bin/env bash
#SBATCH --job-name=answerALS_gvcf_to_txt
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/gvcf_to_txt/wrapper_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/gvcf_to_txt/wrapper_%j.err

set -euo pipefail
export LC_ALL=C

BASE="/home/zw529/donglab/data/answer_ALS"
RNA_ROOT="${BASE}/iPSC_Motor_Neuron/RNAseq/Processed"
WGS_ROOT="${BASE}/Blood/WGS/Processed"
QTL="${BASE}/QTL"
REFERENCE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
OUTPUT_ROOT="${QTL}/chromosome_joint_vcfs/gvcf_to_txt"

module load BCFtools/1.21-GCC-13.3.0

MODE="${1:-submit}"

# Canonical chromosomes only; decoy and alt contigs are excluded.
CHROMOSOMES=( {1..22} X Y M )
N_CHROMOSOMES="${#CHROMOSOMES[@]}"
EXPECTED_DONORS=815

submit_jobs() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "[ERROR] Run submit mode from a login shell with bash, not through sbatch." >&2
        exit 1
    fi

    [[ -s "$REFERENCE" ]] || {
        echo "[ERROR] Missing reference FASTA: $REFERENCE" >&2
        exit 1
    }
    [[ -s "${REFERENCE}.fai" ]] || {
        echo "[ERROR] Missing reference index: ${REFERENCE}.fai" >&2
        exit 1
    }
    [[ -d "$RNA_ROOT" ]] || {
        echo "[ERROR] Missing RNA-seq donor directory: $RNA_ROOT" >&2
        exit 1
    }

    local script run_id run_dir manifest chrom_file sample_map
    local array_response array_id merge_response merge_id final_response final_id
    local rna_dir donor donor_wgs_dir input sample task
    local n_rna_dirs=0 n_missing_gvcf=0
    local -a rna_dirs files source_samples

    script="$(readlink -f "$0")"
    run_id="$(date +%Y%m%d_%H%M%S)_$$"
    run_dir="${OUTPUT_ROOT}/run_${run_id}"

    mkdir -p "$run_dir"/{expanded,joint_by_chr,logs,metrics}
    manifest="$run_dir/donor_manifest.tsv"
    chrom_file="$run_dir/chromosomes.txt"
    sample_map="$run_dir/sample_map.txt"

    printf '%s\n' "${CHROMOSOMES[@]}" > "$chrom_file"
    printf 'task\tdonor\tgvcf_sample_id\tinput_gvcf\n' > "$manifest"
    : > "$sample_map"

    mapfile -t rna_dirs < <(
        find "$RNA_ROOT" -mindepth 1 -maxdepth 1 -type d | sort
    )

    declare -A seen_samples=()
    task=0

    for rna_dir in "${rna_dirs[@]}"; do
        [[ -d "$rna_dir" ]] || continue
        n_rna_dirs=$((n_rna_dirs + 1))
        donor="$(basename "$rna_dir")"
        donor_wgs_dir="$WGS_ROOT/$donor"

        if [[ ! -d "$donor_wgs_dir" ]]; then
            echo "[INFO] $donor: no matching WGS folder; skipping."
            n_missing_gvcf=$((n_missing_gvcf + 1))
            continue
        fi

        shopt -s nullglob
        files=( "$donor_wgs_dir"/*.haplotypeCalls.er.raw.vcf.gz )
        shopt -u nullglob

        if ((${#files[@]} == 0)); then
            echo "[INFO] $donor: no raw gVCF in its WGS folder; skipping."
            n_missing_gvcf=$((n_missing_gvcf + 1))
            continue
        elif ((${#files[@]} != 1)); then
            echo "[ERROR] $donor: expected one raw gVCF; found ${#files[@]}." >&2
            exit 1
        fi

        input="${files[0]}"
        mapfile -t source_samples < <(bcftools query -l "$input")

        if [[ "${#source_samples[@]}" -ne 1 ]]; then
            echo "[ERROR] Expected one sample in $input; found ${#source_samples[@]}." >&2
            exit 1
        fi

        sample="${source_samples[0]}"

        if [[ -n "${seen_samples[$sample]+x}" ]]; then
            echo "[ERROR] Duplicate source sample ID '$sample' for $donor and ${seen_samples[$sample]}." >&2
            exit 1
        fi
        seen_samples["$sample"]="$donor"

        printf '%s\t%s\t%s\t%s\n' \
            "$task" "$donor" "$sample" "$input" >> "$manifest"
        printf '%s\t%s\n' "$sample" "$donor" >> "$sample_map"

        for chr in "${CHROMOSOMES[@]}"; do
            mkdir -p "$run_dir/expanded/chr${chr}"
        done

        task=$((task + 1))
    done

    echo "RNA-seq donor folders examined: $n_rna_dirs"
    echo "RNA-seq donors without a matching raw gVCF: $n_missing_gvcf"
    echo "RNA/WGS donors selected: $task"

    if [[ "$task" -ne "$EXPECTED_DONORS" ]]; then
        echo "[ERROR] Expected $EXPECTED_DONORS donors with both RNA-seq and raw WGS gVCF; found $task." >&2
        echo "[ERROR] Refusing to submit an incomplete or differently sized cohort." >&2
        exit 1
    fi

    echo "Run directory: $run_dir"
    echo "Chromosomes: ${CHROMOSOMES[*]}"
    echo "Expansion array: ${EXPECTED_DONORS}"

    array_response=$(
        sbatch --parsable \
            --time=23:00:00 \
            --cpus-per-task=1 \
            --mem=12G \
            --array="1-$((EXPECTED_DONORS))%95" \
            --job-name=answerALS_expand \
            --output="$run_dir/logs/expand_%A_%a.out" \
            --error="$run_dir/logs/expand_%A_%a.err" \
            "$script" worker "$run_dir"
    )
    array_id="${array_response%%;*}"

    merge_response=$(
        sbatch --parsable \
            --partition=week \
            --time=5-23:00:00 \
            --cpus-per-task=4 \
            --mem=48G \
            --dependency="afterok:${array_id}" \
            --array="1-$((N_CHROMOSOMES))%20" \
            --job-name=answerALS_merge_chr \
            --output="$run_dir/logs/merge_%A_%a.out" \
            --error="$run_dir/logs/merge_%A_%a.err" \
            "$script" merge_chr "$run_dir"
    )
    merge_id="${merge_response%%;*}"

    final_response=$(
        sbatch --parsable \
            --partition=week \
            --time=6-00:00:00 \
            --cpus-per-task=4 \
            --mem=24G \
            --dependency="afterok:${merge_id}" \
            --job-name=answerALS_concat \
            --output="$run_dir/logs/concat_%j.out" \
            --error="$run_dir/logs/concat_%j.err" \
            "$script" concat "$run_dir"
    )
    final_id="${final_response%%;*}"

    echo
    echo "Submitted expansion array: $array_id (815 RNA/WGS donors, %95)"
    echo "Submitted chromosome-merge array: $merge_id (25 chromosomes, %4; after expansion)"
    echo "Submitted final concat/index job: $final_id (after chromosome merges)"
    echo "Outputs will be under: $run_dir"
}

expand_one_donor() {
    local run_dir="$1"
    local task_id="${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
    local manifest="$run_dir/donor_manifest.tsv"
    local row task donor sample input output tmp
    local start records total_records=0 chr

    row="$(awk -F '\t' -v n="$task_id" 'NR == n + 2 {print; exit}' "$manifest")"
    [[ -n "$row" ]] || {
        echo "[ERROR] No manifest row for expansion task $task_id." >&2
        exit 1
    }

    IFS=$'\t' read -r task donor sample input <<< "$row"

    [[ "$task" == "$task_id" ]] || {
        echo "[ERROR] Manifest task mismatch for task $task_id." >&2
        exit 1
    }
    [[ -s "$input" ]] || {
        echo "[ERROR] Missing input gVCF: $input" >&2
        exit 1
    }

    if [[ ! -s "${input}.csi" && ! -s "${input}.tbi" ]]; then
        echo "[$donor] Creating source gVCF index..."
        bcftools index --threads 1 -c "$input"
    fi

    start=$(date +%s)

    for chr in "${CHROMOSOMES[@]}"; do
        output="$run_dir/expanded/chr${chr}/${donor}.chr${chr}.expanded.vcf.gz"
        tmp="${output}.tmp.${SLURM_JOB_ID:-$$}"

        echo "[$donor] Expanding chr${chr}..."

        # Stream expanded VCF text directly into BGZF; do not retain a huge
        # uncompressed text intermediate.
        bcftools convert --gvcf2vcf \
            --fasta-ref "$REFERENCE" \
            --regions "chr${chr}" \
            --output-type v \
            "$input" |
        bgzip -@ 1 -c > "$tmp"

        bcftools index -f -t "$tmp"
        mv "$tmp" "$output"
        mv "${tmp}.tbi" "${output}.tbi"

        [[ "$(bcftools query -l "$output")" == "$sample" ]] || {
            echo "[ERROR] Sample ID changed during expansion for $donor chr${chr}." >&2
            exit 1
        }

        records="$(bcftools index -n "$output")"
        total_records=$((total_records + records))
        echo "[$donor] chr${chr}: records=$records compressed_bytes=$(stat -c %s "$output")"
    done

    printf 'task\tdonor\tsample_id\texpanded_records_all_chromosomes\ttotal_seconds\n%s\t%s\t%s\t%s\t%s\n' \
        "$task" "$donor" "$sample" "$total_records" "$(( $(date +%s) - start ))" \
        > "$run_dir/metrics/task_${task_id}.tsv"

    echo "[$donor] Completed all chromosomes: total_records=$total_records"
}

merge_one_chromosome() {
    local run_dir="$1"
    local chr_id="${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
    local chr manifest="$run_dir/donor_manifest.tsv"
    local sample_map="$run_dir/sample_map.txt"
    local list="$run_dir/lists"
    local merged_dir="$run_dir/joint_by_chr"
    local file old_tmp renamed_tmp output
    local expected actual

    chr="$(sed -n "$((chr_id + 1))p" "$run_dir/chromosomes.txt")"
    [[ -n "$chr" ]] || {
        echo "[ERROR] No chromosome for merge task $chr_id." >&2
        exit 1
    }

    mkdir -p "$list"
    list="$list/chr${chr}.vcfs.list"
    : > "$list"

    while IFS=$'\t' read -r task donor sample input; do
        [[ "$task" == "task" ]] && continue
        file="$run_dir/expanded/chr${chr}/${donor}.chr${chr}.expanded.vcf.gz"
        [[ -s "$file" && ( -s "${file}.tbi" || -s "${file}.csi" ) ]] || {
            echo "[ERROR] Missing expanded VCF or index: $file" >&2
            exit 1
        }
        printf '%s\n' "$file" >> "$list"
    done < "$manifest"

    [[ "$(wc -l < "$list" | tr -d ' ')" -eq "$EXPECTED_DONORS" ]] || {
        echo "[ERROR] Expected $EXPECTED_DONORS inputs for chr${chr}." >&2
        exit 1
    }

    old_tmp="$merged_dir/.chr${chr}.old_ids.vcf.gz"
    renamed_tmp="$merged_dir/.chr${chr}.renamed.vcf.gz"
    output="$merged_dir/answer_ALS_chr${chr}.vcf.gz"

    echo "Coordinate-merging ${EXPECTED_DONORS} expanded gVCFs for chr${chr}..."
    bcftools merge \
        --threads 4 \
        --merge both \
        --file-list "$list" \
        --output-type z \
        --output "$old_tmp"

    bcftools reheader \
        --samples "$sample_map" \
        --output "$renamed_tmp" \
        "$old_tmp"

    bcftools index -f -t "$renamed_tmp"
    mv "$renamed_tmp" "$output"
    mv "${renamed_tmp}.tbi" "${output}.tbi"
    rm -f "$old_tmp"

    expected="$(tail -n +2 "$manifest" | cut -f2 | sort)"
    actual="$(bcftools query -l "$output" | sort)"
    [[ "$expected" == "$actual" ]] || {
        echo "[ERROR] Donor IDs in merged chr${chr} VCF do not match the manifest." >&2
        exit 1
    }

    echo "Completed chr${chr}: samples=$(bcftools query -l "$output" | wc -l | tr -d ' ') records=$(bcftools index -n "$output") bytes=$(stat -c %s "$output")"
}

concatenate_chromosomes() {
    local run_dir="$1"
    local list="$run_dir/joint_chromosomes.list"
    local chr file final_tmp final expected actual

    : > "$list"
    for chr in "${CHROMOSOMES[@]}"; do
        file="$run_dir/joint_by_chr/answer_ALS_chr${chr}.vcf.gz"
        [[ -s "$file" && ( -s "${file}.tbi" || -s "${file}.csi" ) ]] || {
            echo "[ERROR] Missing merged chromosome VCF or index: $file" >&2
            exit 1
        }
        printf '%s\n' "$file" >> "$list"
    done

    final_tmp="$run_dir/.joint_genotyped_GQ.vcf.gz.tmp"
    final="$run_dir/joint_genotyped_GQ.vcf.gz"

    echo "Concatenating chromosome VCFs..."
    bcftools concat \
        --threads 4 \
        --file-list "$list" \
        --output-type z \
        --output "$final_tmp"

    bcftools index -f -t "$final_tmp"
    mv "$final_tmp" "$final"
    mv "${final_tmp}.tbi" "${final}.tbi"

    expected="$(tail -n +2 "$run_dir/donor_manifest.tsv" | cut -f2 | sort)"
    actual="$(bcftools query -l "$final" | sort)"
    [[ "$expected" == "$actual" ]] || {
        echo "[ERROR] Final joint VCF sample IDs do not match the donor manifest." >&2
        exit 1
    }

    echo
    echo "Final sample count: $(bcftools query -l "$final" | wc -l | tr -d ' ')"
    echo "Final record count: $(bcftools index -n "$final")"
    echo "Final compressed bytes: $(stat -c %s "$final")"
    echo "Final joint VCF: $final"
}

case "$MODE" in
    submit)
        submit_jobs
        ;;
    worker)
        [[ $# -eq 2 ]] || { echo "Usage: $0 worker RUN_DIR" >&2; exit 1; }
        expand_one_donor "$2"
        ;;
    merge_chr)
        [[ $# -eq 2 ]] || { echo "Usage: $0 merge_chr RUN_DIR" >&2; exit 1; }
        merge_one_chromosome "$2"
        ;;
    concat)
        [[ $# -eq 2 ]] || { echo "Usage: $0 concat RUN_DIR" >&2; exit 1; }
        concatenate_chromosomes "$2"
        ;;
    *)
        echo "Usage: $0 [submit]" >&2
        exit 1
        ;;
esac
