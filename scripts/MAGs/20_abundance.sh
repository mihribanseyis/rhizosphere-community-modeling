#!/usr/bin/env bash

#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=abund_RhizoMAGs
#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------

if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch $0 GENOME_DIR READS_DIR SAMPLE_LIST OUTPUT_DIR" >&2
    exit 1
fi

genome_dir="$1"
reads_dir="$2"
sample_list="$3"
output_dir="$4"

threads="${SLURM_CPUS_PER_TASK:-8}"

renamed_dir="${output_dir}/renamed_genomes"
combined_reference="${output_dir}/all_RhizoMAGs.fa"
coverage_table="${output_dir}/coverage.tsv"


# ---------------------------------------------------------------------------
# Software environment
# ---------------------------------------------------------------------------

module load bio/BWA/0.7.17-GCC-10.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0

conda_profile="${CONDA_PROFILE:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
coverm_environment="${COVERM_CONDA_ENV:-base}"

if [[ ! -f "$conda_profile" ]]; then
    echo "ERROR: Conda initialization script not found: $conda_profile" >&2
    exit 1
fi

source "$conda_profile"
conda activate "$coverm_environment"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

if [[ ! -d "$genome_dir" ]]; then
    echo "ERROR: Dereplicated MAG directory not found: $genome_dir" >&2
    exit 1
fi

if [[ ! -d "$reads_dir" ]]; then
    echo "ERROR: Filtered-read directory not found: $reads_dir" >&2
    exit 1
fi

if [[ ! -s "$sample_list" ]]; then
    echo "ERROR: Sample list not found or empty: $sample_list" >&2
    exit 1
fi

for program in bwa samtools coverm; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

mapfile -t genome_files < <(
    find "$genome_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.fa' \
        -print \
        | sort
)

mag_count="${#genome_files[@]}"

if [[ "$mag_count" -eq 0 ]]; then
    echo "ERROR: No dereplicated MAGs found in: $genome_dir" >&2
    exit 1
fi

mapfile -t samples < <(
    grep -v '^[[:space:]]*$' "$sample_list" |
        grep -v '^[[:space:]]*#'
)

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "ERROR: No sample IDs found in: $sample_list" >&2
    exit 1
fi

mkdir -p "$output_dir" "$renamed_dir"


# ---------------------------------------------------------------------------
# Rename MAG contig headers
# ---------------------------------------------------------------------------

start_time=$(date +%s)

echo "[$(date)] Starting abundance estimation"
echo "Genome directory: $genome_dir"
echo "MAG count:        $mag_count"
echo "Reads directory:  $reads_dir"
echo "Sample count:     ${#samples[@]}"
echo "Output directory: $output_dir"
echo "Threads:          $threads"
echo

echo "[$(date)] Renaming MAG contig headers"

find "$renamed_dir" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    -delete

for genome_file in "${genome_files[@]}"; do
    mag_id="$(basename "$genome_file" .fa)"
    renamed_file="${renamed_dir}/${mag_id}.fa"

    awk -v prefix="$mag_id" '
        BEGIN {
            OFS = ""
        }

        /^>/ {
            header = substr($0, 2)
            split(header, fields, /[ \t]/)

            remainder = ""

            if (length(header) > length(fields[1])) {
                remainder = substr(header, length(fields[1]) + 1)
            }

            print ">", prefix, "__", fields[1], remainder
            next
        }

        {
            print
        }
    ' "$genome_file" > "$renamed_file"

    if [[ ! -s "$renamed_file" ]]; then
        echo "ERROR: Renamed MAG file was not created or is empty: $renamed_file" >&2
        exit 1
    fi
done

renamed_count=$(
    find "$renamed_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.fa' \
        | wc -l
)

if [[ "$renamed_count" -ne "$mag_count" ]]; then
    echo "ERROR: Expected $mag_count renamed MAGs, but found $renamed_count." >&2
    exit 1
fi

echo "[$(date)] Renamed $renamed_count MAGs"


# ---------------------------------------------------------------------------
# Build combined MAG reference
# ---------------------------------------------------------------------------

echo "[$(date)] Concatenating renamed MAGs"

cat "${renamed_dir}"/*.fa > "$combined_reference"

if [[ ! -s "$combined_reference" ]]; then
    echo "ERROR: Combined MAG reference was not created or is empty." >&2
    exit 1
fi

echo "[$(date)] Building BWA index"

bwa index "$combined_reference"


# ---------------------------------------------------------------------------
# Map filtered reads against the MAG reference
# ---------------------------------------------------------------------------

bam_files=()
missing_samples=0

for sample in "${samples[@]}"; do
    read1="${reads_dir}/${sample}/non_host_R1_${sample}.fastq.gz"
    read2="${reads_dir}/${sample}/non_host_R2_${sample}.fastq.gz"

    bam_file="${output_dir}/${sample}.sorted.bam"
    bam_index="${bam_file}.bai"

    if [[ ! -s "$read1" || ! -s "$read2" ]]; then
        echo "WARNING: Paired reads missing or empty for $sample" >&2
        echo "         R1: $read1" >&2
        echo "         R2: $read2" >&2

        missing_samples=$((missing_samples + 1))
        continue
    fi

    echo "[$(date)] Mapping $sample"

    bwa mem \
        -t "$threads" \
        "$combined_reference" \
        "$read1" \
        "$read2" |
        samtools view \
            -@ "$threads" \
            -b \
            - |
        samtools sort \
            -@ "$threads" \
            -o "$bam_file" \
            -

    samtools index "$bam_file"

    if [[ ! -s "$bam_file" || ! -s "$bam_index" ]]; then
        echo "ERROR: BAM file or index was not created for $sample." >&2
        exit 1
    fi

    bam_files+=("$bam_file")
done

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files were produced." >&2
    exit 1
fi

echo "[$(date)] Produced ${#bam_files[@]} BAM files"


# ---------------------------------------------------------------------------
# Estimate genome coverage with CoverM
# ---------------------------------------------------------------------------

echo "[$(date)] Running CoverM genome"

coverm genome \
    --bam-files "${bam_files[@]}" \
    --genome-fasta-directory "$renamed_dir" \
    --genome-fasta-extension fa \
    --threads "$threads" \
    --output-file "$coverage_table"

if [[ ! -s "$coverage_table" ]]; then
    echo "ERROR: CoverM coverage table was not created or is empty." >&2
    exit 1
fi


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

end_time=$(date +%s)
runtime_minutes=$(( (end_time - start_time) / 60 ))

echo
echo "[$(date)] Abundance estimation completed"
echo "MAGs included:     $mag_count"
echo "Samples requested: ${#samples[@]}"
echo "Samples mapped:    ${#bam_files[@]}"
echo "Samples missing:   $missing_samples"
echo "Coverage table:    $coverage_table"
echo "Combined reference: $combined_reference"
echo "Renamed MAGs:      $renamed_dir"
echo "Runtime:           ${runtime_minutes} minutes"

if [[ "$missing_samples" -gt 0 ]]; then
    echo "WARNING: Some samples were skipped because their reads were missing." >&2
fi