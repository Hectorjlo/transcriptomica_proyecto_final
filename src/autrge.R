#!/usr/bin/env Rscript


extract_attr_value <- function(attributes, key) {
	# Extract the value of one attribute key from the attributes
	# of a gff3 file, if a key is not found returns the original string }
	# and converts them to NA_character_
	# 
	# Returns:
	#	Character vector of extracted values or NA_character_

	# gff3 style: key=value
	gff33_pattern <- paste0(".*(?:^|;)", key, "=([^;]+).*")
	value <- sub(gff33_pattern, "\\1", attributes, perl = TRUE)
	value[value == attributes] <- NA_character_

	value
}

load_tx2gene_from_annotation <- function(path) {
	# Builds a transcript to gene mapping table from a gff3 annotation file.
	# Reads feature type and attributes column of the gff3, keeps rows where the feature
	# type is a transcript or mRNA, extract those transcript IDs and gene IDs, creates
	# a table and cleans it from duplications
	# 
	# Returns:
	#	A two column table of transcript, gene cleaned from duplications
	
	# Reads the gff3 file
	gff3 <- read.delim(
			path,
			sep = "\t",
			header = FALSE,
			stringsAsFactors = FALSE,
			quote = "",
			comment.char = "#",
			colClasses = c(
				"NULL", "NULL", "character", "NULL", "NULL",
				"NULL", "NULL", "NULL", "character"
			)
		)

	# Extracts feature type and attributes column from the file
	feature_type <- gff3[[1]]
	attributes <- gff3[[2]]

	# Select the rows were feature type is transcript or mRNA 
	keep <- feature_type %in% c("transcript", "mRNA")
	attributes <- attributes[keep]

	# Extracts transcript IDs (tx_id)
	tx_id <- extract_attr_value(attributes, "transcript_id") # Searches for transcript_id
	tx_fallback <- extract_attr_value(attributes, "ID") # Fallback to seach ID
	tx_id[is.na(tx_id)] <- tx_fallback[is.na(tx_id)] # Select tx_id or tx_fallback is tx_id is NA
	tx_id <- sub("^transcript:", "", tx_id) # Removes transcript prefix
	tx_id <- sub("\\.[0-9]+$", "", tx_id, perl = TRUE) # Removes the version suffix like .1 or .2 etc

	# Extracts the gene ID
	gene_id <- extract_attr_value(attributes, "gene_id") # Searches for gene_id
	gene_fallback <- extract_attr_value(attributes, "Parent") # Fallback to seach Parent
	gene_id[is.na(gene_id)] <- gene_fallback[is.na(gene_id)] # Select fallback if gene_id is NA
	gene_id <- sub("^gene:", "", gene_id) # Removes prefix

	# Creates a table of two columns, one for the transcript_id and one for the gene
	tx2gene <- unique(data.frame(
		transcript = tx_id,
		gene = gene_id,
		stringsAsFactors = FALSE # Do not convert string to factors
	))

	# Removes rows with empty transcript or gene
	tx2gene <- tx2gene[nzchar(tx2gene$transcript) & nzchar(tx2gene$gene), , drop = FALSE]

	# If a transcript maps to multiple genes, keeps only the first
	if (anyDuplicated(tx2gene$transcript)) {
		warning("Duplicate transcript IDs found in tx2gene mapping. Keeping first occurrence.")
		tx2gene <- tx2gene[!duplicated(tx2gene$transcript), , drop = FALSE]
	}

	# Returns the cleaned table
	tx2gene
}

run_tximport_and_write <- function(quant_files, output_tsv, tx2gene) {
	# Runs transcript import (tximport) and writes the final table
	# Returns:
	#	A table with gene_id, transcript abundance for each quant file


	# Runs tximport on the quant files (salmon generated) and the table tx2gene
	txi <- tximport::tximport(
		files = quant_files,
		type = "salmon",
		tx2gene = tx2gene,
		ignoreTxVersion = TRUE
	)

	# Sets the output dir, creates it if it not found
	out_dir <- dirname(output_tsv)
	if (!dir.exists(out_dir)) {
		dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	}
	# Save the rds tximport object
  rds_path <- sub("\\.tsv$", ".rds", output_tsv)
  saveRDS(txi, file = rds_path)
	# Gathers the gene abundance generated from tx abundance
	gene_abundance <- data.frame(
		gene_id = rownames(txi$abundance),
		txi$abundance,
		check.names = FALSE,
		row.names = NULL,
		stringsAsFactors = FALSE
	)

	# Writes the table to a file
	write.table(
		gene_abundance,
		file = output_tsv,
		sep = "\t",
		quote = FALSE,
		row.names = FALSE
	)

}

# main:

# Parses the arguments of the CLI
args <- commandArgs(trailingOnly = TRUE)
input_dir_single <- args[[1]]
input_dir_paired <- args[[2]]
output_tsv_single <- args[[3]]
output_tsv_paired <- args[[4]]
# Path of the annotation file
annotation_file <- "/home/hectorjl/WORKING_DIR/4to/transcriptomica/DE_analysis/data/gencode.vM36/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz"

# Loads and sort the quant files of single end reads
quant_files_single <- list.files(
	path = input_dir_single,
	pattern = "^quant\\.sf$",
	full.names = TRUE,
	recursive = TRUE
)
quant_files_single <- sort(quant_files_single)

# Gathers by REGEX the core name of the SRR
sample_names_single <- basename(dirname(quant_files_single))
sample_names_single <- sub("^output_", "", sample_names_single)
sample_names_single <- sub("\\.salmon$", "", sample_names_single)
sample_names_single <- make.unique(sample_names_single, sep = "_")
names(quant_files_single) <- sample_names_single

# Load and sort paired end quant files
quant_files_paired <- list.files(
	path = input_dir_paired,
	pattern = "^quant\\.sf$",
	full.names = TRUE,
	recursive = TRUE
)
quant_files_paired <- sort(quant_files_paired)
# Gathers by REGEX core name
sample_names_paired <- basename(dirname(quant_files_paired))
sample_names_paired <- sub("^output_", "", sample_names_paired)
sample_names_paired <- sub("\\.salmon$", "", sample_names_paired)
sample_names_paired <- make.unique(sample_names_paired, sep = "_")
names(quant_files_paired) <- sample_names_paired

# Loads the annotation file
tx2gene <- load_tx2gene_from_annotation(annotation_file)

# Run the function for single end reads
run_tximport_and_write(
	quant_files = quant_files_single,
	output_tsv = output_tsv_single,
	tx2gene = tx2gene
)

# Runs the function for paired end reads
run_tximport_and_write(
	quant_files = quant_files_paired,
	output_tsv = output_tsv_paired,
	tx2gene = tx2gene
)

