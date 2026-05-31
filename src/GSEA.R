# Carga librerías necesarias
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# Carga el objeto RDS generado del análisis diferencial
res <- readRDS(file = "results/star/DESeq2/results.rds")

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
            OrgDb = org.Hs.eg.db, # Base de datos de Humano
            eps = 1e-300 # No trunca el P-value hasta 1e-300
        )

# Guarda el objeto RDS generado
saveRDS(gse, file = "results/star/GSEA/gse.rds")

# Realiza el plot
p <- gseaplot2(gse,
               geneSetID = c("GO:0030198", "GO:1901342"),
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
ggsave(filename = "results/star/GSEA/plots/GSEA_plot.png",
    plot = p, 
    dpi = 1200, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)
