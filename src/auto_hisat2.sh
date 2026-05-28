#!/usr/bin/env bash

# Runs hisat2 single_end in parallel for each read, harcoded code.
# This script expects two arguments, one that is the path in which the
# files will be search on and the second where output files will be generated
# Usage:
#   ./auto_hisat2.sh <input_dir> <output_dir>
#   ./auto_hisat2.sh ../data/trimmed_fastqs ../results/hisat2/single_end

# Same principle as `auto_fastp.sh`
# Create an array of the names that match: clean_*_1.fastq
mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "clean_*_1.fastq" | sort)

# Count how many read_files are
n_files=$(( ${#files[@]} ))


# Hardcoded function
run_hisat2() {
    # Arguments
    #   -x: Path of the index of the reference genome 
    #   -U: Input of the fastq file of single_end reads, taken from the first argument passed to the function ($1)
    #   -S: Path of the output SAM file ($2) 
    #   --summary-file: Path of output summary file ($3)
    #   -k: Maximum number of alignments reported per read 
    #   --no-unal: Don't report unaligned reads in the SAM file
    hisat2 -x /export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat \
           -U "$1" \
           -S "$2" \
           --summary-file "$3" \
           -k 7 \
           --no-unal

}
# Export the function, needed for parallel to use it
export -f run_hisat2

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
#   ../../../data/trimmed_fastqs/cleaned_SRR1111.fastq --> SRR1111
core_name="${files[$i]##*/}"
core_name="${core_name#clean_}"
core_name="${core_name%.fastq}"
# The echo will send 3 strings:
#   1st -> The fastp file of unpaired reads 
#   2nd -> Name of the output SAM file
#   3rd -> Name of the summary file 
# These strings are build in the echo
echo "${files[$i]} ${out_dir}output_${core_name}.sam ${out_dir}summary_of_${core_name}.txt"
done | parallel -j ${n_files} --colsep ' ' run_hisat2