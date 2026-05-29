#!/usr/bin/env bash

# Runs featureCounts on all BAM files found in the input directory simultaneously.
# This produces a single unified count matrix and avoids reloading the GFF3 file multiple times.
#
# Usage:
#   ./auto_ft_counts.sh <input_dir> <output_file>
# Example:
#   ./auto_ft_counts.sh ../results/bam ../results/counts/all_counts.txt

INPUT_DIR="$1"
OUTPUT_FILE="$2"
ANNOTATION="/home/hectorjl/4to/transcriptomica_proyecto_final/data/GENCODE_GRCh38.p13_104/gencode.v38.chr_patch_hapl_scaff.annotation.gff3"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Find all BAM files in input directory
mapfile -t bam_files < <(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.bam" | sort)

# Run featureCounts on all BAM files at once using 24 threads
featureCounts -o "$OUTPUT_FILE" \
    -a "$ANNOTATION" \
    -T 24 \
    --largestOverlap \
    "${bam_files[@]}"

