# Comentario
## Documentación

###
# Uso: Rscript DEA_DESeq2.R <gene_counts_or_txi_file> <annotacion_file> <gene_map_file> <result_dir> 
###

# Desactiva el dispositivo PDF
## Evita generar archivos fuera de los indicados
pdf(file = NULL)

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
    })
)

# Genera un parseador
parser <- list(
    # Obtiene el path de counts o txi
    make_option(
        "--counts",
        type = "character"
    ),
    # Obtiene el path del archivo de anotación
    make_option(
        "--annotation",
        type = "character"
    ),
    # Obtiene el path del archivo gene_map (Geneid - gene_name) 
    make_option(
        "--gene_map",
        type = "character"
    ),
    # Obtiene el path de la carpeta resultados
    make_option(
        "--results_dir",
        type = "character"
    ),
    # Bandera para diferenciar si se procesara un objeto txi
    # o una tabla de conteos
    make_option(
        "--from_pseudoalignment",
        action = "store_true",
        default = FALSE
    )
)

# Parsea el objeto para ser accedido por índice
args <- parse_args(OptionParser(option_list = parser))

# Obten los argumentos del paseo de la CLI
gene_counts_or_txi_path <- args$counts
annotacion_file_path <- args$annotation
gene_name_map_file_path <- args$gene_map
results_files_dir <- args$results_dir
from_pseudoalignment <- args$from_pseudoalignment

## Para el path de "results" elimina si encuentra a un "/"
## y agrega un "/", asegurando una única aparición 
results_files_dir <- gsub("/*$", "/", results_files_dir)

# Lee los archivos de anotación y del genemap
annotation <- read.delim(annotacion_file_path, row.names=1)
gene_name_map <- read.delim(gene_name_map_file_path, header=FALSE, row.names=1)

# Secuencia ordenada de cabeceras de los datos en procesamiento
#!
reordered <- c("male_24m_1", "male_24m_2", "male_24m_4", "male_24m_7",
                "male_3m_3", "male_3m_5", "male_3m_6")

# Si es un objeto txi
if (from_pseudoalignment) {
    # Lee el objeto txi
    gene_counts_or_txi <- readRDS(gene_counts_or_txi_path)
    # Registra el nombre de cada SRR respecto a la condición
    #!
    condition <- c("male_24m_2", "male_24m_1", "male_3m_5", "male_24m_7", "male_3m_6", "male_3m_3", "male_24m_4")
    # Renombra los nombres de las columnas
    colnames(gene_counts_or_txi$abundance) <- condition
    colnames(gene_counts_or_txi$counts) <- condition
    colnames(gene_counts_or_txi$length) <- condition
    # Reordena las columnas por edad
    gene_counts_or_txi$abundance <- gene_counts_or_txi$abundance[ ,reordered]
    gene_counts_or_txi$counts <- gene_counts_or_txi$counts[ ,reordered]
    gene_counts_or_txi$length <- gene_counts_or_txi$length[ ,reordered]
} else { # Si es una tabla de conteos
    # Lee el archivo de conteos
    gene_counts_or_txi <- read.delim(gene_counts_or_txi_path, row.names=1)
    # Reordena las columnas por edad
    gene_counts_or_txi <- gene_counts_or_txi[ ,reordered]
}

# Generar la tabla de metadatos
# Verifica con colnames(gene_counts_or_txi), comprueba la cabecera de los archivos
## Para estos archivos es:
## "male_24m_1" "male_24m_2" "male_24m_4" "male_24m_7" "male_3m_3"  "male_3m_5"  "male_3m_6"
## Por cual, primero crea un factor con el mismo orden que la cabecera en "age"
age <- factor(c(rep("m24", 4), rep("m3", 3)),
    levels = c("m3", "m24"))
## Al igual que otro por "sex" siguiendo la cabecera de "sex"
sex <- factor(c(rep("male", 7)))
## Además, un factor que por color refleje a cual de los cuatro grupos pertencen
## lightblue <- male_3m
## blue <- male_24m
sample_color <- c(rep("blue", 4), rep("lightblue", 3))
# Si es un objeto txi
if (from_pseudoalignment) {
    ## Por último los colnames() 
    sample_names <- colnames(gene_counts_or_txi$counts)
} else { # Si es una tabla de conteos
    ## Por último los colnames() 
    sample_names <- colnames(gene_counts_or_txi)
}
# Con los valores anteriores genera la tabla de metadatos
meta_data <- data.frame(sample_names, age, sex)
## Genera
##   sample_names age  sex
## 1   male_24m_1 m24 male
## 2   male_24m_2 m24 male
## 3   male_24m_4 m24 male
## 4   male_24m_7 m24 male
## 5    male_3m_3  m3 male
## 6    male_3m_5  m3 male
## 7    male_3m_6  m3 male
# Elimina la columna "sample_names", pero usa su contenido
# como nombre de filas
meta_data <- meta_data %>% 
        remove_rownames %>% 
        column_to_rownames(var="sample_names") 
## Genera
##            age  sex
## male_24m_1 m24 male
## male_24m_2 m24 male
## male_24m_4 m24 male
## male_24m_7 m24 male
## male_3m_3   m3 male
## male_3m_5   m3 male
## male_3m_6   m3 male

# Verifica que en los nombres de columnas en gene_counts_or_txi
# esten el mismo orden que los nombred de filas en meta_data
print("Verificación: colnames(gene_counts_or_txi) == rownames(meta_data)")
print(ifelse(
    (all(colnames(gene_counts_or_txi) == rownames(meta_data))),
    "Ok",
    "Error"
))

# Crea un objeto de DESeq2 usando los datos anteriores
# Usa de diseño de matriz a: ~ 0 + age
## Usando "~ 0 + age", se usa un modelo de medias de grupo
## lo cual permite poder crear los contrastes de comparación
## de manera más fina
formula <- ~ 0 + age

# Si es un objeto txi
if (from_pseudoalignment) {
    # Utiliza la formula para objetos de tximport
    dds <- DESeqDataSetFromTximport(
        txi = gene_counts_or_txi, 
        colData = meta_data, 
        design = formula
    )
} else { # Utiliza la formula para una tabla de conteos
    dds <- DESeqDataSetFromMatrix(
        countData = round(gene_counts_or_txi),
        colData = meta_data,
        design = formula
    )
}

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
# Convierte a vst
vsd <- vst(dds)
# Genera el PCA plot
PCA_plot <- plotPCA(vsd, intgroup = "age") +  
    theme_minimal(base_size = 18, base_line_size = 1)

# Guarda el gráfico de PCA
ggsave(filename = paste(results_files_dir, "plots/PCA_plot.png", sep = ""), 
    plot = PCA_plot, 
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
contrasts <- makeContrasts(m24_vs_m3 = agem24 - agem3,
    levels = design
)

# Realiza el análisis de expresión diferencial
## En este caso al ser solo uno es directo
## Se usa el contraste comparando agem24 vs agem3
res <- results(dds, contrast = contrasts[ ,"m24_vs_m3"])
# Añade Gene name
res$Gene_name <- gene_name_map[rownames(res), ]

# Fija los umbrales de "Log Fold Change" y "False Discovery Rate"
FDR <- 5e-2 # 0.05 FDR
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
    row.names = T
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
    row.names = T
)

# Forma el "volcano plot"
# Asigna colores
volcano_plot_colors <- c("gray", "blue", "red")
names(volcano_plot_colors) <- c("NO", "DOWN", "UP")

# Añade la columna DE
res$DE <- "NO"
res[up, "DE"] <- "UP"
res[down, "DE"] <- "DOWN"

volcano_plot <- ggplot(
    data = res,
    aes(x = log2FoldChange, y = -log10(padj), col = DE)
) +
    geom_point() +
    labs(
        title = "Volcano plot",
        x = "log2 Expression fold change", 
        y = "-log10 FDR"
    ) + 
    scale_color_manual(values = volcano_plot_colors) +
    geom_vline(
        xintercept = c(-LogFC, LogFC),
        col = "black",
        linetype = "longdash"
    ) +
    geom_hline(
        yintercept = -log10(FDR),
        col = "black",
        linetype = "longdash"
    ) +
    theme_minimal(
        base_size = 18,
        base_line_size = 1
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

# Reordena las muestras por edad
age_order <- order(meta_data$age)
# Reordena z-score según la edad
z_score_significant_ordered <- z_score_significant[, age_order]

# Plotea
heatmap_plot <- Heatmap(
    z_score_significant,
    cluster_rows = TRUE,
    cluster_columns = FALSE, 
    show_row_names = FALSE,
    name = "Z score",
    km = 2,
    column_title = "Heatmap all",
    column_labels = colnames(z_score_significant_ordered),
    col = colorRamp2(c(-2, -1, 0, 1, 2), 
        c("#07f900", "#007500", "#000000", "#750b00", "#ff2500"))
)
# Guarda el mapa de calor
png(paste(results_files_dir, "plots/heatmap_all.png", sep = ""), 
    width = 11.25, height = 7.5, res = 300, units = "in")
draw(heatmap_plot)
dev.off()

# Formar Heatmap top20
significant <- head(rownames(significant_order), n = 20)

# Guardar la tabla de los top 20 genes significativos
write.table(
    significant_order[significant, ],
    file = paste(results_files_dir, "DEG_top20_table.tsv", sep = ""),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)

z_score_top_20 <- t(scale(t(log2_tmp[significant, ])))
top_20_heatmap <- Heatmap(
    z_score_top_20,
    cluster_rows = T, 
    cluster_columns = F, 
    row_labels = gene_name_map[rownames(z_score_top_20), ], 
    name = "Z-score", 
    km = 2, 
    column_title = "Top 20 significant genes",
    col = colorRamp2(c(-2, -1, 0, 1, 2), 
        c("#07f900", "#007500", "#000000", "#750b00", "#ff2500"))
)
# Guarda el mapa de calor
png(paste(results_files_dir, "plots/heatmap_top20.png", sep = ""), 
    width = 11.25, height = 7.5, res = 300, units = "in")
draw(top_20_heatmap)
dev.off()

# Guarda el RDS de res para analisis tipo GSEA/RNAAgeCalc
saveRDS(res, file = paste(results_files_dir, "results.rds", sep = ""))

# Cierra dispositivos gráficos
graphics.off()