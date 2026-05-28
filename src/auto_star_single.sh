#!/usr/bin/env bash

# Runs STAR single_end in parallel for each read, hardcoded code.
# This script expects two arguments, one that is the path in which the
# files will be search on and the second where output files will be generated
# Usage:
#   ./auto_star_single.sh <input_dir> <output_dir>
#   ./auto_star_single.sh ../data/trimmed_fastqs ../results/star/single_end

# Same principle as `auto_fastp.sh` and `auto_hisat2.sh` and `auto_hisat2_pe.sh`
# Create an array of the names that match: clean_*_1.fastq
mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*_1.fastq" | sort)

# Count how many read_files are
n_files=$(( ${#files[@]} ))

# Hardcoded function
run_star() {
    # Arguments
    #   --runMode alignReads: Set STAR to alignment mode
    #   --genomeDir: Path of the STAR index of the reference genome
    #   --readFilesIn: Input of the fastq file of single_end reads, taken from the first argument passed to the function ($1)
    #   --outSAMtype SAM: Output alignment file in SAM format
    #   --outFileNamePrefix: Prefix for STAR output files ($2)
    #   --outSAMunmapped None: Don't include unmapped reads in SAM output, similar to --no-unal in hisat2
    #   --runThreadN: Number of threads for STAR execution
    #   --outFilterMultimapNmax: Maximum number of multiple alignments allowed per read, similar to -k in hisat2
    STAR --runMode alignReads \
         --genomeDir /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.star/ \
         --readFilesIn "$1" \
         --outSAMtype SAM \
         --outFileNamePrefix "$2" \
         --outSAMunmapped None \
         --runThreadN 2 \
         --outFilterMultimapNmax 7
}
# Export the function, needed for parallel to use it
export -f run_star

# Take the second argument passed to the main call ($2) and remove a "/" adding at the end a "/"
out_dir="${2%/}/"
# For each read_file do ... until i is less than the number of elements in the array "files"
# starting with i = 0
for (( i=0; i<${#files[@]}; i+=1 )); do
# Build the core_name var, it takes the element "i" of the array "files"
# starting from the front of the string will eliminate all "*" until a "/" is found,
# then again starting from the front will eliminate the first instance of "clean_",
# finally starting from the back of the string will eliminate the ".fastq"
# e.g:
#   ../../../data/trimmed_fastqs/clean_SRR1111_1.fastq --> SRR1111_1
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%.fastq}"
# The echo will send 2 strings:
#   1st -> The fastq file of single-end reads
#   2nd -> STAR output prefix path
# These strings are build in the echo
echo "${files[$i]} ${out_dir}output_${core_name}_"
done | parallel -j ${n_files} --colsep ' ' run_star