#!/usr/bin/env bash

# Runs STAR paired_end in parallel for each read, hardcoded code.
# This script expects two arguments, one that is the path in which the
# files will be search on and the second where output files will be generated
# Usage:
#   ./auto_star_pe.sh <input_dir> <output_dir>
#   ./auto_star_pe.sh ../data/trimmed_fastqs ../results/star/paired_end

# Same principle as `auto_fastp.sh` and `auto_hisat2.sh` and `auto_hisat2_pe.sh` and `auto_star_single.sh`
# Create an array of the names that match: clean_*.fastq
mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*.fastq" | sort)

# Count how many read_files are / 2, pe reads
n_files=$(( ${#files[@]} / 2 ))

# Hardcoded function
run_star() {
    # Arguments
    #   --runMode alignReads: Set STAR to alignment mode
    #   --genomeDir: Path of the STAR index of the reference genome
    #   --readFilesIn: Input of paired-end reads, taken from the first and second arguments passed to the function ($1 and $2)
    #   --outSAMtype SAM: Output alignment file in SAM format
    #   --outFileNamePrefix: Prefix for STAR output files ($3)
    #   --outSAMunmapped None: Don't include unmapped reads in SAM output, similar to --no-unal in hisat2
    #   --runThreadN: Number of threads for STAR execution
    #   --outFilterMultimapNmax: Maximum number of multiple alignments allowed per read, similar to -k in hisat2
    #   --outSAMattributes NH HI AS NM: Include key alignment tags in the SAM output
    STAR --runMode alignReads \
         --genomeDir /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.star/ \
         --readFilesIn "$1" "$2" \
         --outSAMtype SAM \
         --outFileNamePrefix "$3" \
         --outSAMunmapped None \
         --runThreadN 2 \
         --outFilterMultimapNmax 7 \
         --outSAMattributes NH HI AS NM 
}
# Export the function, needed for parallel to use it
export -f run_star

# Take the second argument passed to the main call ($2) and remove a "/" adding at the end a "/"
out_dir="${2%/}/"
# For each read_file do ... until i is less than the number of elements in the array "files"
# starting with i = 0, advancing by 2 because each sample has _1 and _2
for (( i=0; i<${#files[@]}; i+=2 )); do
# Build the core_name var, it takes the element "i" of the array "files"
# starting from the front of the string will eliminate all "*" until a "/" is found,
# then again starting from the front will eliminate the first instance of "clean_",
# finally starting from the back of the string will eliminate the "_1.fastq"
# e.g:
#   ../../../data/trimmed_fastqs/clean_SRR1111_1.fastq --> SRR1111
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%_1.fastq}"
# The echo will send 3 strings:
#   1st -> The fastq file of forward reads
#   2nd -> The fastq file of reverse reads
#   3rd -> STAR output prefix path
# These strings are build in the echo
echo "${files[$i]} ${files[$i+1]} ${out_dir}output_${core_name}_pe_"
done | parallel -j ${n_files} --colsep ' ' run_star