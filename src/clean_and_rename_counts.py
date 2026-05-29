#!/usr/bin/env python3
import os
import re

# Paths
samples_file = "data/GSE213001/extracted_samples.tsv"
input_counts = "results/star/feature_counts/ft_counts_matrix.tsv"
output_counts = "results/star/feature_counts/counts_matrix.tsv"

def main():
    # 1. Read the samples metadata and map Run -> Diagnosis
    run_to_diagnosis = {}

    with open(samples_file, "r") as f:
        header = f.readline().strip().split("\t")
        run_idx = header.index("Run")
        diag_idx = header.index("Diagnosis")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > max(run_idx, diag_idx):
                run_id = parts[run_idx]
                diagnosis = parts[diag_idx]
                run_to_diagnosis[run_id] = diagnosis

    print(f"Loaded {len(run_to_diagnosis)} samples from metadata.")

    with open(input_counts, "r") as infile, open(output_counts, "w") as outfile:
        # Skip the comment line (first line starting with #)
        first_line = infile.readline()
        if not first_line.startswith("#"):
            # If it wasn't a comment line, we need to process it as the header
            header_line = first_line
        else:
            header_line = infile.readline()

        header_cols = header_line.strip().split("\t")
        
        # We need the 1st column (Geneid) and columns from index 6 (7th column) onwards
        bam_mappings = []

        for idx, col in enumerate(header_cols):
            if idx >= 6:
                # Extract SRR identifier (e.g. from results/star/alignments/output_SRR21498140_Aligned.out.bam)
                match = re.search(r"SRR\d+", col)
                if match:
                    run_id = match.group(0)
                    diagnosis = run_to_diagnosis.get(run_id, "Unknown")
                    new_name = f"{diagnosis}_{run_id}"
                    bam_mappings.append((diagnosis, run_id, idx, new_name))
                else:
                    print(f"Warning: Could not extract SRR ID from column: {col}")
                    bam_mappings.append(("Unknown", col, idx, col))

        # Sort the mappings: IPF first, then NDC (alphabetical works because 'I' < 'N'), and sub-sort by run_id
        bam_mappings.sort(key=lambda x: (x[0] != "IPF", x[0] != "NDC", x[0], x[1]))

        # Build the sorted new header and get corresponding sorted indices
        new_header = [header_cols[0]] + [item[3] for item in bam_mappings]
        bam_cols_indices = [item[2] for item in bam_mappings]

        # Write the new header to output
        outfile.write("\t".join(new_header) + "\n")

        # Process each data row
        count = 0
        for line in infile:
            cols = line.strip().split("\t")
            if len(cols) > 0:
                gene_id = cols[0]
                row_data = [gene_id]
                for idx in bam_cols_indices:
                    row_data.append(cols[idx])
                outfile.write("\t".join(row_data) + "\n")
                count += 1

    print(f"Successfully processed {count} genes.")
    print(f"Output saved to: {output_counts}")

if __name__ == "__main__":
    main()
