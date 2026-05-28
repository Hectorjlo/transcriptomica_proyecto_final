#!/usr/bin/env bash
# Before running a check CPU usage, this creates a job using 2 threads
# per file all at the same time so if 7 files are found it would use
# 7*2 = 14 threads
# Executes a auto processing of SAM files to ordered BAM files
# This script expects two arguments, one that is the path in which the
# files will be search on and the second where output files will be generated
# Usage:
#   ./auto_BAMic.sh <input_dir> <output_dir>
#   ./auto_BAMic.sh ../results/bam ../results/sam

# Create an array of the names that match: *.sam
mapfile -t files < <(find "$1" -maxdepth 1 -type f -name "*.sam" | sort)
# Calculate total file length
n_files=$(( ${#files[@]} ))

# Hardcoded function
run_samtools() {
    # Arguments
    #   -@: Threads per job
    #   -b: Output in BAM format
    #   -S: Input type file auto detect
    #   -q: Keep alignments with MAPQ >= 10
    #   -T: Temporal files names (autodeleted when the processing ends)
    #   -o: Output BAM file path
    samtools view -@ 2 -bSq 10 "$1"| samtools sort -@ 2 -T "$2" -o "$3"

}
# Export the function, needed for parallel to use it
export -f run_samtools

# Take the second argument passed to the main call ($2) and remove a "/" adding at the end a "/"
out_dir="${2%/}/"
# For each read_file do ... until i is less than the number of elements in the array "files"
# starting with i = 0
for (( i=0; i<${#files[@]}; i+=1 )); do
# Build the core_name var, it takes the element "i" of the array "files"
# starting from the front of the string will eliminate all "*" until a "/" is found,
# then again starting from the front will eliminate all ultil a "SRR" match is found,
# finally starting from the back of the string will eliminate the ".sam"
# e.g:
#   ../../../results/output_SRR9126859_pe_Aligned.out.sam --> SRR9126859_pe_Aligned.out
core_name="${files[$i]##*/}"
core_name="SRR${core_name#*SRR}"
core_name="${core_name%.sam}"
# The echo will send 3 strings:
#   1st -> The input file
#   2nd -> The core name
#   3rd -> The outfile path
# These strings are build in the echo
echo "${files[$i]} ${core_name} ${out_dir}output_${core_name}.bam"
done | parallel -j ${n_files} --colsep ' ' run_samtools