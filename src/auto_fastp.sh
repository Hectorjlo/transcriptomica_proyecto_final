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
    #   -i/-I: Read 1 and Read 2 of pair end files ($1/$2)
    #   -o/-O: Trimmed output files sufix, prefix is "clean_" ($3/$4)
    #   --adapter_sequence/--adapter_sequence_r2: An adapter secuence to be found 
    #                                             and deleted when files are processed
    #   --trim_front1/--trim_front2: Number of bases to be trimmed in the front of each read
    #   --detect_adapter_for_pe: Detect and delete common adapters
    #   --trim_poly_g: Trim sequence errors (Common in the platform used for this fastq)
    #   --l: Minimum length for read after trimming
    fastp -i "$1" \
          -o "$2" \
          --trim_front1 12 \
          -D \
          --dup_calc_accuracy 5 \
          --trim_poly_g \
          --trim_poly_x -l 50
}
# Export a function (-f) to the environment, enables the call to the function
export -f run_fastp

# Adds an "/" in case the path does not contain one, and keeps only one in the other case
# ../output/dir -> ../output/dir/
# /path/outdir/ -> /path/outdir/
out_dir="${2%/}/"

# For each do ... until i is less than the number of elements in the array "files"
# starting with i = 0
for (( i=0; i<${#files[@]}; i+=1 )); do
#                           ^^^ Pass to the next 
# The echo send 4 strings, the first two are the path of the files fastq
# the next two are only the names of those files (deletes the */ part keeps the name)
    echo "${files[$i]} ${out_dir}${files[$i]##*/}" 
done | parallel -j "$n_files" --colsep ' ' run_fastp  
#      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ For each (separated by spaces) execute the function
# By this way each fastp is run in a separate theread, concurrency at it most


