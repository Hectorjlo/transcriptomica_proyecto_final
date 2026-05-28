# Carga librerías necesarias
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Carga el objeto RDS generado del análisis diferencial
res <- readRDS(file = "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/results.rds")

# Ordena de manera inversa por stat
res <- res[order(-res$stat), ]

# Crea la lista de genes
gene_list <- res$stat
# Renombra el vector igual que los nombres de los genes del objeto RDS
names(gene_list) <- rownames(res) 

# Realiza el análisis de GSE
gse <- gseGO(geneList = gene_list,
            ont = "BP", # Agrega la ontología "Biological Process"
            keyType = "ENSEMBL", # Formato ENSEMBL ID
            OrgDb = "org.Mm.eg.db", # Base de datos de ratón
            eps = 1e-300 # No trunca el P-value hasta 1e-300
        )

# Guarda el objeto RDS generado
saveRDS(gse, file = "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/gse.rds")

# Realiza el plot
p <- gseaplot2(gse,
               geneSetID = c(grep(
                "osteoclast proliferation", 
                as.data.frame(gse)$Description
                ), # Busca a ese patrón retorna el índice
                grep("negative regulation of osteoclast differentiation", 
                as.data.frame(gse)$Description)
                ),
               title = "",
               base_size = 14,
               pvalue_table = TRUE, # Muestra la tabla de P-value
               pvalue_table_columns = c("NES", "p.adjust") # Agrega el valor NES
            )

p[[1]] <- p[[1]] +  geom_hline( # Agrega una línea en y = 0
    yintercept = 0, 
    linetype = "solid", 
    color = "grey40", 
    linewidth = 0.5
)

# Guarda el plot
ggsave(filename = "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/plots/GSEA_plot.png",
    plot = p, 
    dpi = 1200, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)
