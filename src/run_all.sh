#!/usr/bin/env bash

### DESeq2
# hisat2
Rscript src/DEA_DESeq2.R --counts results/hisat2/feature_counts/single_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/hisat2/DE_analysis/DESeq2/single_end
Rscript src/DEA_DESeq2.R --counts results/hisat2/feature_counts/paired_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/hisat2/DE_analysis/DESeq2/paired_end

# star
Rscript src/DEA_DESeq2.R --counts results/star/feature_counts/single_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/star/DE_analysis/DESeq2/single_end
Rscript src/DEA_DESeq2.R --counts results/star/feature_counts/paired_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/star/DE_analysis/DESeq2/paired_end
# salmon
Rscript src/DEA_DESeq2.R --counts results/salmon/tximport/single_end/gene_ab.rds \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/salmon/DE_analysis/DESeq2/single_end --from_pseudoalignment
Rscript src/DEA_DESeq2.R --counts results/salmon/tximport/paired_end/gene_ab.rds \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/salmon/DE_analysis/DESeq2/paired_end --from_pseudoalignment

### EdgeR
Rscript src/DEA_EdgeR.R --counts results/hisat2/feature_counts/single_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/hisat2/DE_analysis/edgeR/single_end
Rscript src/DEA_EdgeR.R --counts results/hisat2/feature_counts/paired_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/hisat2/DE_analysis/edgeR/paired_end

# star
Rscript src/DEA_EdgeR.R --counts results/star/feature_counts/single_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/star/DE_analysis/edgeR/single_end
Rscript src/DEA_EdgeR.R --counts results/star/feature_counts/paired_end/feature_counts_global.tsv \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/star/DE_analysis/edgeR/paired_end
# salmon
Rscript src/DEA_EdgeR.R --counts results/salmon/tximport/single_end/gene_ab.rds \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/salmon/DE_analysis/edgeR/single_end --from_pseudoalignment
Rscript src/DEA_EdgeR.R --counts results/salmon/tximport/paired_end/gene_ab.rds \
        --annotation data/gencode.vM36/gencode.vM36.gene_id.length.tsv --gene_map data/gencode.vM36/mm39-gencode-M36-gene_id-gene_name.txt \
        --results_dir results/salmon/DE_analysis/edgeR/paired_end --from_pseudoalignment