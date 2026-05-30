# Comentario
## Documentación

###
# Uso: Rscript DEA_DESeq2.R <matrix_counts_file> <annotacion_file> <gene_map_file> <result_dir> 
###


# Carga de las librerías necesarias para el análisis
# Suprime Warnings para una salida stdout más limpia
suppressWarnings(
    suppressPackageStartupMessages({
        library(DESeq2)
        library(ggplot2)
        library(ComplexHeatmap)
        library(dplyr)
        library(tibble)
        library(edgeR)
        library(circlize)
        library(optparse)
        library(ggrepel)
    })
)

# Carga los archivos
matrix_counts_path <- "results/star/feature_counts/counts_matrix.tsv"
annotacion_file_path <- "results/star/feature_counts/gene_id.gene_length.tsv"
gene_name_map_file_path <- "data/GENCODE_GRCh38.p13_104/gene_id.gene_name.txt"
results_files_dir <- "results/star/DESeq2"
from_pseudoalignment <- FALSE
samples_table_path <- "data/GSE213001/extracted_samples.tsv"


## Para el path de "results" elimina si encuentra a un "/"
## y agrega un "/", asegurando una única aparición 
results_files_dir <- gsub("/*$", "/", results_files_dir)

# Lee los archivos de anotación y del genemap
annotation <- read.delim(annotacion_file_path, row.names=1)
gene_name_map <- read.delim(gene_name_map_file_path, header=FALSE, row.names=1)
# Elimina versión del geneid de ENSEMBL
rownames(gene_name_map) <- gsub("[.][0-9]+","",rownames(gene_name_map))

# Se cargí la tabla de muestras
samples_table <- read.delim(samples_table_path)
## Debido a que un dato fué encontrado como NA, pero el PCA de edad
## mostró que no había un efecto de lote, se decidio tomar la media
## para asignarle una edad y mantener los datos que ofrece
samples_table$AGE[is.na(samples_table$AGE)] <- 61
# Carga la matriz de conteos
matrix_counts <- read.delim(matrix_counts_path, row.names=1)
# Elimina versión del ENSEMBL geneid 
rownames(matrix_counts) <- gsub("[.][0-9]+","",rownames(matrix_counts))

# Elimina la ruta, mantiene el ID de SRR
run_ids <- sub(".*(SRR[0-9]+).*", "\\1", colnames(matrix_counts))
# Ordena la tabla de muestras dado el ID de la matriz de conteos
idx <- match(run_ids, samples_table$Run)

# Cambia los nombres para reflejar su origen
colnames(matrix_counts) <- paste(
  samples_table$diagnosis[idx],
  samples_table$gender[idx],
  samples_table$AGE[idx],
  run_ids,
  sep = "_"
)

# Ordena por condición
diagnosis_order <- factor(
  samples_table$diagnosis[idx],
  levels = c("IPF", "NDC")
)

# Ordena la matriz de conteos dado la condición
matrix_counts <- matrix_counts[, order(diagnosis_order)]

# Generar la tabla de metadatos
## Dado el ordenamiento anterior fue 20 IPF y 14 NDC
condition <- factor(c(rep("IPF", 20), rep("NDC", 14)),
    levels = c("NDC", "IPF"))

# Crea el factor así como se encuentre "male" o "female" en las columnas 
sex <- factor(stringr::str_extract(colnames(matrix_counts), "(?<=_)(male|female)(?=_)"),
    levels = c("male", "female"))

# Crea el factor de age así como se encuentra en las columnas
age <- stringr::str_extract(colnames(matrix_counts), "(?<=_)(\\d+|NA)(?=_)") |>
  as.numeric()

# Agrega color para distinguir IPF de NDC
sample_color <- c(rep("blue", 12), rep("lightblue", 12))

sample_names <- colnames(matrix_counts)
# Con los valores anteriores genera la tabla de metadatos
meta_data <- data.frame(sample_names, condition, sex, age)

# Elimina la columna "sample_names", pero usa su contenido
# como nombre de filas
meta_data <- meta_data %>% 
        remove_rownames %>% 
        column_to_rownames(var="sample_names") 

# Verifica que en los nombres de columnas en matrix_counts
# esten el mismo orden que los nombred de filas en meta_data
print("Verificación: colnames(matrix_counts) == rownames(meta_data)")
print(ifelse(
    (all(colnames(matrix_counts) == rownames(meta_data))),
    "Ok",
    "Error"
))

# Crea un objeto de DESeq2 usando los datos anteriores
# Usa de diseño de matriz a: ~ 0 + condition
## Usando "~ 0 + condition", se usa un modelo de medias de grupo
## lo cual permite poder crear los contrastes de comparación
## de manera más fina
formula <- ~ 0 + condition + sex + age

# Se crea el objeto de DESeq
dds <- DESeqDataSetFromMatrix(
    countData = round(matrix_counts),
    colData = meta_data,
    design = formula
)

# Crea la matriz de diseño
design <- model.matrix(formula)

# Muestra genes totales
print("Genes totales:")
print(length(rownames(dds)))

# Filtrar los genes con baja expresión
keep <- filterByExpr(dds, design)
print("Genes después de filtrado por expresión:")
print(sum(keep))
# Selecciona solo los genes que pasan el filtro
dds <- dds[keep, ]
# Elimina el filtro por uso de memoria
rm(keep)

## Se creará el PCA plot, para ello se debe 
## convertir los datos a vst
vsd <- vst(dds)
# Genera el PCA plot de condición
PCA_condition <- plotPCA(vsd, intgroup = "condition") +  
    theme_minimal(base_size = 18, base_line_size = 1)

# PCA por edad
PCA_age <- plotPCA(vsd, intgroup = "age") +
  theme_minimal(base_size = 18, base_line_size = 1)

# PCA por sexo
PCA_sex <- plotPCA(vsd, intgroup = "sex") +
  theme_minimal(base_size = 18, base_line_size = 1)



# Guarda el gráfico de PCA
ggsave(filename = paste(results_files_dir, "plots/PCA_condition.png", sep = ""), 
    plot = PCA_condition, 
    dpi = 300, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)

ggsave(filename = paste(results_files_dir, "plots/PCA_age.png", sep = ""), 
    plot = PCA_age, 
    dpi = 300, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)

ggsave(filename = paste(results_files_dir, "plots/PCA_sex.png", sep = ""), 
    plot = PCA_sex, 
    dpi = 300, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)

# Calcula los factores de normalización,
# varianza y ajusta el modelo
dds <- DESeq(dds)

## Se calculó los TPM, métrica importante pero no indespensable
# Añade la longitud de los genes
mcols(dds)$basepairs <- annotation[rownames(dds), ]
# Calcula el log_2 FPKM
## Se añadió un pseudo conteo para evitar log2(0)
log2_fpkm <- log2(fpkm(dds) + 0.1)
# Función, convierte de fpkm a log2(tmp)
fpkm_to_tpm_log2 <- function(fpkm) {
    fpkm - log2(sum(2^fpkm)) + log2(1e6)
} 
# Aplica la función por columnas (2)
log2_tmp <- apply(log2_fpkm, 2, fpkm_to_tpm_log2)
# Guarda en una tabla
gene_names <- gene_name_map[rownames(log2_tmp), ]
write.table(
    cbind(
        gene_names,
        log2_tmp
    ),
    file = paste(results_files_dir, "TPM_log2-table.txt", sep = ""),
    sep = "\t",
    quote = FALSE
)

## Se realizó los contrastes de expresión diferencial
# Imprime nombres disponibles para los contrastes
print("Nombres disponibles para contrastes:")
print(resultsNames(dds))

# Crea el contraste
contrasts <- makeContrasts(IPF_vs_NDC = conditionIPF - conditionNDC,
    levels = design
)

# Realiza el análisis de expresión diferencial
## En este caso al ser solo uno es directo
## Se usa el contraste comparando conditionIPF vs conditionNDC
res <- results(dds, contrast = contrasts[ ,"IPF_vs_NDC"])
# Añade Gene name
res$Gene_name <- gene_name_map[rownames(res), ]

# Fija los umbrales de "Log Fold Change" y "False Discovery Rate"
FDR <- 0.001
LogFC <- 0.5

# Obten los genes "upregulated" y "downregulated"
up <- (res$log2FoldChange > LogFC) & (res$padj < FDR)
# Elimina los NA
up[which(is.na(up))] <- FALSE
# Imprimir los "upregulated"
cat("Upregulated: ", sum(up), "\n")

down <- (res$log2FoldChange < -LogFC) & (res$padj < FDR)
# Elimuna los NA
down[which(is.na(down))] <- FALSE
cat("Downregulated: ", sum(down), "\n")

# Guarda el resumen de genes diferenciales en un archivo
writeLines(
    c(
        paste("Upregulated:", sum(up)),
        paste("Downregulated:", sum(down))
    ),
    con = paste(
            results_files_dir,
            "DEG_summary.txt",
            sep = ""
        )
)

# Guardar las tablas para "up" y "down" regulated
write.table(
    res[up, ],
    paste(results_files_dir, "deseq-DEG-up", 
        FDR, 
        ".txt", 
        sep =""
    ),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)
write.table(
    res[down, ],
    paste(results_files_dir, "deseq-DEG-down", 
        FDR, 
        ".txt", 
        sep = ""
    ),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)

# Forma el "volcano plot"
# Asigna colores
volcano_plot_colors <- c("#d9d9d9", "#e06666", "#93c47d")
names(volcano_plot_colors) <- c("NO", "DOWN", "UP")

# Añade la columna DE
res$DE <- "NO"
res[up, "DE"] <- "UP"
res[down, "DE"] <- "DOWN"
res$label <- ifelse(res$padj < 1e-6 & abs(res$log2FoldChange) > 1, res$Gene_name, "")

volcano_plot <- ggplot(
    data = res,
    aes(x = log2FoldChange, y = -log10(padj), col = DE)
) +
    geom_point(alpha = 0.6, size = 2) +
    labs(
        title = "IPF contra NDC",
        x = "log2 Expression fold change", 
        y = "-log10 FDR"
    ) + 
    scale_color_manual(values = volcano_plot_colors) +
    scale_fill_manual(values = volcano_plot_colors) +
    geom_vline(
        xintercept = c(-LogFC, LogFC),
        col = "black",
        linetype = "longdash"
    ) +
    geom_hline(
        yintercept = -log10(FDR),
        col = "darkgray",
        linetype = "dashed",
        linewidth = 0.5
    ) +
    theme_light(
        base_size = 14,
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "plain", size = 16), 
        panel.background = element_rect(fill = "#f8f9fa"), 
        legend.position = "bottom",
        legend.title = element_text(face = "bold")
    ) + 
    geom_label_repel(
        aes(label = label, fill = DE),
        color = "black",
        segment.color = "grey50",
        segment.size = 0.3,
        box.padding = 0.4,
        point.padding = 0.3,
        max.overlaps = 30,
        show.legend = FALSE,
        alpha = 0.5
    )
# Guarda el gráfico de volcán
ggsave(filename = paste(results_files_dir, "plots/vulcano_plot.png", sep = ""), 
    plot = volcano_plot, 
    dpi = 300, 
    width = 11.25, 
    height = 7.5,
    create.dir = TRUE
)

# Formar el Heatmap
# Obten genes "up" o "down" *regulated*
significant <- res[up | down, ]
# Ordena por FDR
significant_order <- significant[order(significant$padj), ]
# Toma los primeros 2000
significant <- head(rownames(significant_order), n = 2000)
# Calcula el z-score
z_score_significant <- t(scale(t(log2_tmp[significant, ])))

# Reordena las muestras por condición
condition_order <- order(meta_data$condition)
# Reordena z-score según la condición
z_score_significant_ordered <- z_score_significant[, condition_order]

# Plotea
heatmap_plot <- Heatmap(
    z_score_significant_ordered,
    cluster_rows = TRUE,
    cluster_columns = FALSE, 
    show_row_names = FALSE,
    name = "Z score",
    km = 2,
    column_title = "Heatmap all",
    column_labels = colnames(z_score_significant_ordered),
    col = colorRamp2(c(-2, -1, 0, 1, 2), 
        c("#ff2500", "#750b00", "#000000", "#007500", "#07f900"))
)
# Guarda el mapa de calor
png(paste(results_files_dir, "plots/heatmap_all.png", sep = ""), 
    width = 11.25, height = 7.5, res = 300, units = "in")
draw(heatmap_plot)
dev.off()

# Formar Heatmap top30
significant <- head(rownames(significant_order), n = 30)

# Guardar la tabla de los top 30 genes significativos
write.table(
    significant_order[significant, ],
    file = paste(results_files_dir, "DEG_top30_table.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)

z_score_top_30 <- t(scale(t(log2_tmp[significant, ])))
z_score_top_30_ordered <- z_score_top_30[, condition_order]
top_30_heatmap <- Heatmap(
    z_score_top_30_ordered,
    cluster_rows = TRUE, 
    cluster_columns = FALSE, 
    row_labels = gene_name_map[rownames(z_score_top_30_ordered), ], 
    column_labels = colnames(z_score_top_30_ordered),
    name = "Z-score", 
    km = 2, 
    column_title = "Top 30 significant genes",
    col = colorRamp2(c(-2, -1, 0, 1, 2), 
        c("#ff2500", "#750b00", "#000000", "#007500", "#07f900"))
)
# Guarda el mapa de calor
png(paste(results_files_dir, "plots/heatmap_top30.png", sep = ""), 
    width = 11.25, height = 7.5, res = 300, units = "in")
draw(top_30_heatmap)
dev.off()

# Guarda el RDS de res para analisis tipo GSEA/RNAconditionCalc
saveRDS(res, file = paste(results_files_dir, "results.rds", sep = ""))

# Cierra dispositivos gráficos
graphics.off()
