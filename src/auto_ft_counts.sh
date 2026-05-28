#!/usr/bin/env bash

# Runs bamCoverage for paired_end files of expected star and hisat2 alignments
# This script expects eight arguments, 4 inputs, 4 outputs
# Usage:
#   ./auto_ft_counts.sh <input_star_single_path> \
#                       <input_star_paired_path> \
#                       <input_hisat_single_path> \
#                       <input_hisat_paired_path> \
#                       <output_star_single_path> \
#                       <output_star_paired_path> \
#                       <output_hisat_single_path> \
#                       <output_hisat_paired_path> \
#   ./auto_ft_counts.sh ../data/trimmed_fastqs ../results/star/paired_end

# Allow extension of glob pattern
shopt -s extglob

# Create arrays of the names that match: *.bam files
mapfile -t star_single < <(find "$1" -maxdepth 1 -type f -name "*.bam" | sort) # For Star_single
mapfile -t star_paired < <(find "$2" -maxdepth 1 -type f -name "*.bam" | sort) # For Star_paired
mapfile -t hisat_single < <(find "$3" -maxdepth 1 -type f -name "*.bam" | sort) # For Hisat2_single
mapfile -t hisat_paired < <(find "$4" -maxdepth 1 -type f -name "*.bam" | sort) # For Hisat2_paired
all_arrays=(star_single star_paired hisat_single hisat_paired)
# Count how many files where found
n_files_star=$(( ${#star_single[@]} )) # For Star
n_files_hisat=$(( ${#hisat_single[@]} )) # For Hisat2


run_ft_counts() {
    # Arguments
    #   -o "$2": Name of output file including read counts
    #   -a: Name of an annotation file. GTF/GFF format
    #   -T: Number of the threads
    #   -a: bin size of 20 bp
    #   --largestOverlap: Assign reads to a meta-feature/feature that has the largest number of overlapping bases
    #   -p: Pair ends reads
    #   -B: Only count read pairs that have both ends aligned
    #   "$1": Path of the BAM file
    if $3; then
        featureCounts -o "$2" \
         -a /home/hectorjl/4to/transcriptomica/data/gencode.vM36/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz \
         -T 4 \
         --largestOverlap \
         "$1"    
    else
    featureCounts -o "$2" \
         -a /home/hectorjl/4to/transcriptomica/data/gencode.vM36/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz \
         -T 4 \
         --largestOverlap \
         -p \
         -B \
         "$1"
    fi
}
# Export the function, needed for parallel to use it
export -f run_ft_counts


# Gather file output paths
out_dir_star_se="${5%/}/" # For Star_se
out_dir_star_pe="${6%/}/" # For Star_pe
out_dir_hisat_se="${7%/}/" # For Hisat_se
out_dir_hisat_pe="${8%/}/" # For Hisat_pe


process_files() {
    # Arguments
    #   local -n input_array: Array of BAM files
    #   local output_dir: Output path
    #   local jobs: Jobs executed at the same time
    local -n input_array=$1
    local output_dir=$2
    local jobs=$3

    for ((i = 0; i < ${#input_array[@]}; i++)); do
        core_name="${input_array[$i]##*/}"
        core_name="${core_name#output_}"
        core_name="${core_name%_pe@(_Aligned.out.bam|.bam)}"

        echo "${input_array[$i]} ${output_dir}ft_counts_${core_name}.txt $4"
    done | parallel -j "$jobs" --colsep ' ' run_ft_counts
}

process_files star_single  "$out_dir_star_se"   "$n_files_star" true
process_files star_paired  "$out_dir_star_pe"   "$n_files_star" false
process_files hisat_single "$out_dir_hisat_se"  "$n_files_hisat" true
process_files hisat_paired  "$out_dir_hisat_pe"  "$n_files_hisat" false

