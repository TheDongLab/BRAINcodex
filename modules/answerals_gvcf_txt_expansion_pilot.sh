#!/bin/bash
#SBATCH --job-name=als_gvcf_expand_pilot
#SBATCH --output=/home/zw529/donglab/data/answer_ALS/QTL/gvcf_expand_chr22_pilot_%j.out
#SBATCH --error=/home/zw529/donglab/data/answer_ALS/QTL/gvcf_expand_chr22_pilot_%j.err
#SBATCH --partition=week
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=32G

set -euo pipefail
export LC_ALL=C

module load BCFtools/1.21-GCC-13.3.0

ROOT="/home/zw529/donglab/data/answer_ALS/Blood/WGS/Processed"
QTL="/home/zw529/donglab/data/answer_ALS/QTL"
REFERENCE="/home/zw529/donglab/references/genome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

RUN_ID="${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}"
PILOT_DIR="$QTL/chromosome_joint_vcfs/expanded_chr22_pilot_${RUN_ID}"
mkdir -p "$PILOT_DIR"/{expanded,logs}

[[ -s "$REFERENCE" ]] || { echo "[ERROR] Missing reference: $REFERENCE" >&2; exit 1; }
[[ -s "${REFERENCE}.fai" ]] || { echo "[ERROR] Missing reference index: ${REFERENCE}.fai" >&2; exit 1; }

# Randomly select five donor folders that contain a raw gVCF.
mapfile -t DONOR_DIRS < <(
    find "$ROOT" -mindepth 2 -maxdepth 2 -type f \
        -name '*.haplotypeCalls.er.raw.vcf.gz' -printf '%h\n' |
    sort -u | shuf -n 5
)
[[ "${#DONOR_DIRS[@]}" -eq 5 ]] || {
    echo "[ERROR] Found only ${#DONOR_DIRS[@]} eligible donor folders; need 5." >&2
    exit 1
}

LIST="$PILOT_DIR/expanded_vcfs.list"
MANIFEST="$PILOT_DIR/donor_manifest.tsv"
METRICS="$PILOT_DIR/expansion_metrics.tsv"
printf 'donor_folder\tgvcf_sample_id\tinput_gvcf\n' > "$MANIFEST"
printf 'donor_folder\texpanded_records\ttext_vcf_bytes\tcompressed_vcf_bytes\texpand_and_compress_seconds\n' > "$METRICS"
: > "$LIST"

# Resolve inputs and save the exact five-donor selection before launching workers.
declare -a INPUTS=()
declare -a DONORS=()
declare -a SAMPLES=()
declare -a OUTPUTS=()

for donor_dir in "${DONOR_DIRS[@]}"; do
    donor=$(basename "$donor_dir")
    files=( "$donor_dir"/*.haplotypeCalls.er.raw.vcf.gz )
    [[ "${#files[@]}" -eq 1 ]] || {
        echo "[ERROR] Expected exactly one raw gVCF for $donor; found ${#files[@]}." >&2
        exit 1
    }

    input="${files[0]}"
    mapfile -t source_samples < <(bcftools query -l "$input")
    [[ "${#source_samples[@]}" -eq 1 ]] || {
        echo "[ERROR] Expected one sample in $input; found ${#source_samples[@]}." >&2
        exit 1
    }

    sample="${source_samples[0]}"
    output="$PILOT_DIR/expanded/${donor}.chr22.expanded.vcf.gz"

    INPUTS+=( "$input" )
    DONORS+=( "$donor" )
    SAMPLES+=( "$sample" )
    OUTPUTS+=( "$output" )

    printf '%s\t%s\t%s\n' "$donor" "$sample" "$input" >> "$MANIFEST"
    printf '%s\n' "$output" >> "$LIST"
done

echo "Pilot directory: $PILOT_DIR"
echo "Selected donors and source sample IDs:"
cat "$MANIFEST"
echo "Filesystem space before expansion:"
df -h "$PILOT_DIR"

expand_one() {
    local donor="$1"
    local sample="$2"
    local input="$3"
    local output="$4"
    local text_vcf="${output%.gz}"
    local start end text_bytes compressed_bytes records

    start=$(date +%s)

    # Region access needs an index; create one only if neither index exists.
    if [[ ! -s "${input}.csi" && ! -s "${input}.tbi" ]]; then
        bcftools index --threads 1 -c "$input"
    fi

    echo "[$donor] Expanding chr22 from $sample ..."
    bcftools convert --gvcf2vcf \
        --fasta-ref "$REFERENCE" \
        --regions chr22 \
        --output-type v \
        --output "$text_vcf" \
        "$input"

    text_bytes=$(stat -c %s "$text_vcf")

    # Compress the text VCF for indexing and coordinate-aware merging.
    bgzip -@ 1 -f "$text_vcf"
    bcftools index -f -t "$output"

    end=$(date +%s)
    compressed_bytes=$(stat -c %s "$output")
    records=$(bcftools index -n "$output")

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$donor" "$records" "$text_bytes" "$compressed_bytes" "$((end - start))" \
        > "$PILOT_DIR/logs/${donor}.metrics.tsv"

    echo "[$donor] records=$records text_bytes=$text_bytes compressed_bytes=$compressed_bytes seconds=$((end - start))"
}

# Expand the five donors concurrently, one CPU per donor.
PIDS=()
for i in "${!DONORS[@]}"; do
    expand_one "${DONORS[$i]}" "${SAMPLES[$i]}" "${INPUTS[$i]}" "${OUTPUTS[$i]}" &
    PIDS+=( "$!" )
done

failed=0
for pid in "${PIDS[@]}"; do
    if ! wait "$pid"; then
        failed=1
    fi
done
[[ "$failed" -eq 0 ]] || {
    echo "[ERROR] At least one donor expansion failed; stopping before merge." >&2
    exit 1
}

for donor in "${DONORS[@]}"; do
    cat "$PILOT_DIR/logs/${donor}.metrics.tsv" >> "$METRICS"
done

MERGED="$PILOT_DIR/chr22.expanded_five_donor_coordinate_merged.vcf.gz"
echo "Merging the five expanded files by genomic position..."
bcftools merge \
    --threads 4 \
    --merge both \
    --file-list "$LIST" \
    --output-type z \
    --output "$MERGED"

bcftools index -f -t "$MERGED"

echo
echo "=== Pilot results ==="
cat "$METRICS"
printf 'Merged samples: '
bcftools query -l "$MERGED" | wc -l
printf 'Merged chr22 records: '
bcftools index -n "$MERGED"
printf 'Merged VCF compressed bytes: '
stat -c %s "$MERGED"
echo "Pilot files and metrics: $PILOT_DIR"
