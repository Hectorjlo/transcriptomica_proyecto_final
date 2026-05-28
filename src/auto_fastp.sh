#!/usr/bin/env bash

# Runs fastp in parallel for n-fastq files with hardcoded options
# This script expects two arguments, one that is the path in which the
# files will be search on and the second where output files will be generated
# Usage: (fastp needs to be callable)
#   ./auto_fastp.sh <input_dir> <output_dir>
#   ./auto_fastp.sh ../data/SRRs/ ../data/trimmed_fastqs/ 

# Using mapfile to create an array of the files
# find and sort are used to keep only files but folders
mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "*.fastq" | sort)
#          ^^^^ Name of the array

# To run parallel it should take n of the fastq files
n_files=$(( ${#files[@]} ))
#             ^^^^^^^^^ Accessing the number of elements in the array
            
# Function to execute fastp for a pair
run_fastp() {
    # Arguments:
    #   -i: Read file ($1)
    #   -o: Output file name ($2)
    #   --trim_front1: Number of bases to be trimmed in the front of each read
    #   --trim_poly_g: Trim g sequence errors
    #   --trim_poly_x: Trim x sequence errors
    #   --l: Minimum length for read after trimming
    fastp -i "$1" \
          -o "$2" \
          --trim_front1 12 \
          --trim_poly_g \
          --trim_poly_x -l 50
}
# Export a function (-f) to the environment, enables the call to the function
export -f run_fastp

# Adds an "/" in case the path does not contain one, and keeps only one in the other case
out_dir="${2%/}/"

for (( i=0; i<${#files[@]}; i+=1 )); do
# The echo send 2 strings
# input_file
# output_file
    echo "${files[$i]} ${out_dir}${files[$i]##*/}" 
done | parallel -j "$n_files" --colsep ' ' run_fastp  
