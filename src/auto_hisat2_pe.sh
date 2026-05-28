#!/usr/bin/env bash

# Runs hisat2 paired_end in parallel for each read, harcoded code.
# This script expects two arguments, one that is the path in which the
# files will be search on and the second where output files will be generated
# Usage:
#   ./auto_hisat2_pe.sh <input_dir> <output_dir>
#   ./auto_hisat2_pe.sh ../data/trimmed_fastqs ../results/hisat2/paired_end

# Same principle as `auto_hisat2.sh`
# Create an array of the names that match: clean_*.fastq
mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*.fastq" | sort)

# Count how many read_files are / 2, pe reads
n_files=$(( ${#files[@]} / 2 ))


# Hardcoded function
run_hisat2() {
    # Arguments
    #   -x: Path of the index of the reference genome
    #   -1: Input of forward reads, taken from the first argument passed to the function ($1)
    #   -2: Input of reverse reads, taken from the second argument passed to the function ($2)
    #   -S: Path of the output SAM file ($3)
    #   --summary-file: Path of output summary file ($4)
    #   -k: Maximum number of alignments reported per read
    #   --no-unal: Don't report unaligned reads in the SAM file
    #   --no-mixed: Don't report mixed alignments for paired reads
    #   --no-discordant: Don't report discordant alignments for paired reads
    hisat2 -x /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat \
           -1 "$1" \
           -2 "$2" \
           -S "$3" \
           --summary-file "$4" \
           -k 7 \
           --no-unal \
           --no-mixed \
           --no-discordant

}
# Export the function, needed for parallel to use it
export -f run_hisat2

# Take the second argument passed to the main call ($2) and remove a "/" adding at the end a "/"
# For each read_file do ... until i is less than the number of elements in the array "files"
# starting with i = 0
out_dir="${2%/}/"
for (( i=0; i<${#files[@]}; i+=2 )); do
# Build the core_name var, it takes the element "i" of the array "files"
# starting from the front of the string will eliminate all "*" until a "/" is found,
# then again starting from the front will eliminate the first instance of "clean_",
# finally starting from the back of the string will eliminate the "_1.fastq"
# e.g:
#   ../../../data/trimmed_fastqs/clean_SRR1111_1.fastq --> SRR1111
# The echo will send 4 strings:
#   1st -> The fastq file of forward reads
#   2nd -> The fastq file of reverse reads
#   3rd -> Name of the output SAM file
#   4th -> Name of the summary file
# These strings are build in the echo
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%_1.fastq}"
echo "${files[$i]} ${files[$i+1]} ${out_dir}output_${core_name}_pe.sam ${out_dir}summary_of_${core_name}_pe.txt"
done | parallel -j ${n_files} --colsep ' ' run_hisat2