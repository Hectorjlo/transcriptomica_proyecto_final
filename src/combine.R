# Se realizó el parseo de los argumentos
args <- commandArgs(trailingOnly = TRUE)
dir_hisat_se <- args[[1]]
dir_hisat_pe <- args[[2]]
dir_star_se <- args[[3]]
dir_star_pe <- args[[4]]


# Se usó el "glob pattern" para encontrar todos los archivos
paths_hisat_se <- Sys.glob(paste(dir_hisat_se, "ft_counts_SRR*.txt", sep = ""))
paths_hisat_pe <- Sys.glob(paste(dir_hisat_pe, "ft_counts_SRR*.txt", sep = ""))
paths_star_se <- Sys.glob(paste(dir_star_se, "ft_counts_SRR*.txt", sep = ""))
paths_star_pe <- Sys.glob(paste(dir_star_pe, "ft_counts_SRR*.txt", sep = ""))

create_combined_table <- function(paths, outdir) {
    # Se inicializó un par de listas
    list_col7 <- list()
    names <- list()
    
    # Se obtuvo en "table" la primera columna del primer archivo
    table <- read.table(paths[1], sep = "\t")
    # Se retiró el punto y numero después de la identificación del ID
    table <- gsub("\\.[0-9]+", "", table$V1[-1])
    
  # Se iteró por cada archivo de paths
    for (i in 1:length(paths)) {
        # Se leyó y se obtuvo la septima columna
        df <- read.table(file = paths[i], sep = "\t")$V7
        # Se sustituyó para solo obtener el ID del SRR
        df[1] <- gsub(".bam", "", gsub(".*SRR", "SRR", df[1]))
        # Se guardó en "name" el SRR
        name <- df[1]
        # Se guardó en la lista el nombre
        names <- append(names, name)
        
        # Se obtuvo el df empezando desde la segunda fila
        col7 <- df[-1]
        # Se añadió a la lista la columna obtenida
        list_col7 <- append(list_col7, list(col7))
    }
    
    # Se transformó a un data.frame la lista de columnas
    final_table <- as.data.frame(list_col7)
    # Se añadió la columna de Geneid a la tabla final
    final_table <- cbind(table = table, final_table)
    # Se renombraron las columnas para reflejar el Geneid y los SRR
    colnames(final_table) <- c("Geneid", names)
    
    # Se guardó la tabla en el mismo path de donde provienen los archivos
    write.table(
        final_table,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        file = paste(outdir, "feature_counts_global.tsv", sep = "")
    )
}

# Se llamó a la función para cada estado
create_combined_table(paths_hisat_se, dir_hisat_se)
create_combined_table(paths_hisat_pe, dir_hisat_pe)
create_combined_table(paths_star_se, dir_star_se)
create_combined_table(paths_star_pe, dir_star_pe)

print("Hecho.")