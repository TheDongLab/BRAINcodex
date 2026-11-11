#!/usr/bin/env bash
#SBATCH --job-name=als_gvcf_chr22_pilot
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/gvcf_chr22_pilot_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/gvcf_chr22_pilot_%j.err

set -euo pipefail
export LC_ALL=C

ROOT="/home/zw529/donglab/data/answer_ALS/Blood/WGS/Processed"
QTL="/home/zw529/donglab/data/answer_ALS/QTL"
REFERENCE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

module load BCFtools/1.21-GCC-13.3.0

MODE="${1:-submit}"

submit_jobs() {
    # Run this mode from a login shell with: bash prepare_joint_vcf_expand_chr22_pilot.sh
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "[ERROR] Submit mode must be run from a login shell, not inside a Slurm job." >&2
        exit 1
    fi

    [[ -s "$REFERENCE" ]] || { echo "[ERROR] Missing reference: $REFERENCE" >&2; exit 1; }
    [[ -s "${REFERENCE}.fai" ]] || { echo "[ERROR] Missing reference index: ${REFERENCE}.fai" >&2; exit 1; }

    local script run_id run_dir manifest vcf_list
    script="$(readlink -f "$0")"
    run_id="$(date +%Y%m%d_%H%M%S)_$$"
    run_dir="$QTL/chromosome_joint_vcfs/expanded_chr22_pilot_${run_id}"

    mkdir -p "$run_dir"/{expanded,logs,metrics}

    manifest="$run_dir/donor_manifest.tsv"
    vcf_list="$run_dir/expanded_vcfs.list"
    printf 'task\tdonor\tgvcf_sample_id\tinput_gvcf\texpanded_vcf\n' > "$manifest"
    : > "$vcf_list"

    mapfile -t donor_dirs < <(
        find "$ROOT" -mindepth 2 -maxdepth 2 -type f \
            -name '*.haplotypeCalls.er.raw.vcf.gz' -printf '%h\n' |
        sort -u | shuf -n 5
    )

    [[ "${#donor_dirs[@]}" -eq 5 ]] || {
        echo "[ERROR] Found ${#donor_dirs[@]} eligible donors; need 5." >&2
        exit 1
    }

    declare -A seen_samples=()
    local task donor input sample output
    task=0

    for donor_dir in "${donor_dirs[@]}"; do
        donor="$(basename "$donor_dir")"
        shopt -s nullglob
        local files=( "$donor_dir"/*.haplotypeCalls.er.raw.vcf.gz )
        shopt -u nullglob

        [[ "${#files[@]}" -eq 1 ]] || {
            echo "[ERROR] Expected exactly one raw gVCF for $donor; found ${#files[@]}." >&2
            exit 1
        }

        input="${files[0]}"
        mapfile -t source_samples < <(bcftools query -l "$input")
        [[ "${#source_samples[@]}" -eq 1 ]] || {
            echo "[ERROR] Expected one VCF sample in $input; found ${#source_samples[@]}." >&2
            exit 1
        }

        sample="${source_samples[0]}"
        [[ -z "${seen_samples[$sample]+x}" ]] || {
            echo "[ERROR] Duplicate gVCF sample ID among selected donors: $sample" >&2
            exit 1
        }
        seen_samples["$sample"]=1

        output="$run_dir/expanded/${donor}.chr22.expanded.vcf.gz"
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "$task" "$donor" "$sample" "$input" "$output" >> "$manifest"
        printf '%s\n' "$output" >> "$vcf_list"
        task=$((task + 1))
    done

    echo "Pilot output directory: $run_dir"
    echo "Five selected donors:"
    cat "$manifest"

local array_args merge_args

array_args=(
    --parsable
    --time=23:55:00
    --cpus-per-task=1
    --mem=12G
    --array=0-4%5
    --job-name=ALS_chr22_expand
    --output="$run_dir/logs/expand_%A_%a.out"
    --error="$run_dir/logs/expand_%A_%a.err"
)
array_response=$(sbatch "${array_args[@]}" "$script" worker "$run_dir")
array_id="${array_response%%;*}"

merge_args=(
    --parsable
    --time=23:55:00
    --cpus-per-task=2
    --mem=24G
    --dependency="afterok:${array_id}"
    --job-name=ALS_chr22_merge5
    --output="$run_dir/logs/merge_%j.out"
    --error="$run_dir/logs/merge_%j.err"
)
merge_response=$(sbatch "${merge_args[@]}" "$script" merge "$run_dir")
merge_id="${merge_response%%;*}"

    echo
    echo "Submitted five-task array: $array_id"
    echo "Submitted dependent merge job: $merge_id"
    echo "The merge job starts only after all five array tasks succeed."
    echo "Pilot directory: $run_dir"
}

expand_one_donor() {
    local run_dir="$1"
    local task_id="${SLURM_ARRAY_TASK_ID:?Missing SLURM_ARRAY_TASK_ID}"
    local manifest="$run_dir/donor_manifest.tsv"
    local row task donor sample input output text_vcf
    local start expand_seconds text_bytes compress_start compress_seconds total_seconds
    local compressed_bytes records

    row="$(awk -F '\t' -v n="$task_id" 'NR == n + 2 {print; exit}' "$manifest")"
    [[ -n "$row" ]] || { echo "[ERROR] No manifest row for array task $task_id." >&2; exit 1; }

    IFS=$'\t' read -r task donor sample input output <<< "$row"
    [[ "$task" == "$task_id" ]] || { echo "[ERROR] Manifest task mismatch." >&2; exit 1; }
    [[ -s "$input" ]] || { echo "[ERROR] Missing input gVCF: $input" >&2; exit 1; }

    if [[ ! -s "${input}.csi" && ! -s "${input}.tbi" ]]; then
        echo "[$donor] Creating source gVCF index..."
        bcftools index --threads 1 -c "$input"
    fi

    text_vcf="${output%.gz}"
    start=$(date +%s)

    echo "[$donor] Expanding chr22 to per-position text VCF..."
    bcftools convert --gvcf2vcf \
        --fasta-ref "$REFERENCE" \
        --regions chr22 \
        --output-type v \
        --output "$text_vcf" \
        "$input"

    expand_seconds=$(( $(date +%s) - start ))
    text_bytes=$(stat -c %s "$text_vcf")

    compress_start=$(date +%s)
    bgzip -@ 1 -f "$text_vcf"
    compress_seconds=$(( $(date +%s) - compress_start ))

    bcftools index -f -t "$output"
    records=$(bcftools index -n "$output")
    compressed_bytes=$(stat -c %s "$output")
    total_seconds=$(( $(date +%s) - start ))

    [[ "$(bcftools query -l "$output")" == "$sample" ]] || {
        echo "[ERROR] Sample ID changed during expansion for $donor." >&2
        exit 1
    }

    printf 'donor\texpanded_records\ttext_vcf_bytes\tcompressed_vcf_bytes\texpand_seconds\tcompress_seconds\ttotal_seconds\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$donor" "$records" "$text_bytes" "$compressed_bytes" \
        "$expand_seconds" "$compress_seconds" "$total_seconds" \
        > "$run_dir/metrics/task_${task_id}.tsv"

    echo "[$donor] records=$records text_bytes=$text_bytes compressed_bytes=$compressed_bytes expand_seconds=$expand_seconds compress_seconds=$compress_seconds"
}

merge_five_donors() {
    local run_dir="$1"
    local vcf_list="$run_dir/expanded_vcfs.list"
    local manifest="$run_dir/donor_manifest.tsv"
    local merged="$run_dir/chr22.expanded_five_donor_coordinate_merged.vcf.gz"

    [[ -s "$vcf_list" ]] || { echo "[ERROR] Missing expanded VCF list: $vcf_list" >&2; exit 1; }

    while IFS= read -r vcf; do
        [[ -s "$vcf" && ( -s "${vcf}.tbi" || -s "${vcf}.csi" ) ]] || {
            echo "[ERROR] Expanded VCF or index missing: $vcf" >&2
            exit 1
        }
    done < "$vcf_list"

    echo "Coordinate-merging the five expanded chr22 VCFs..."
    bcftools merge \
        --threads 2 \
        --merge both \
        --file-list "$vcf_list" \
        --output-type z \
        --output "$merged"

    bcftools index -f -t "$merged"

    local expected actual samples records
    expected="$(tail -n +2 "$manifest" | cut -f3 | sort)"
    actual="$(bcftools query -l "$merged" | sort)"
    [[ "$expected" == "$actual" ]] || {
        echo "[ERROR] Merged sample IDs do not match the five selected gVCFs." >&2
        exit 1
    }

    samples="$(bcftools query -l "$merged" | wc -l | tr -d ' ')"
    records="$(bcftools index -n "$merged")"

    echo
    echo "=== Expansion measurements ==="
    for metrics_file in "$run_dir"/metrics/task_*.tsv; do
        [[ -s "$metrics_file" ]] && cat "$metrics_file"
    done

    echo
    echo "Merged sample count: $samples"
    echo "Merged chr22 record count: $records"
    echo "Merged compressed bytes: $(stat -c %s "$merged")"
    echo "Merged pilot VCF: $merged"
}

case "$MODE" in
    submit)
        submit_jobs
        ;;
    worker)
        [[ $# -eq 2 ]] || { echo "Usage: $0 worker RUN_DIR" >&2; exit 1; }
        expand_one_donor "$2"
        ;;
    merge)
        [[ $# -eq 2 ]] || { echo "Usage: $0 merge RUN_DIR" >&2; exit 1; }
        merge_five_donors "$2"
        ;;
    *)
        echo "Usage: $0 [submit]" >&2
        exit 1
        ;;
esac
