# Carga las librerías necesarias
library(ggplot2)
library(dplyr)
library(readr)

# Carga los datos de DAVID
david_data <- read_delim(
    "DE_analysis/data/DAVID_files/DAVIDChartReport_upregulated_CC.csv",
    delim = ","
)

# Limpiar y preparar los datos
# Vamos a ordenar por P-Value y quedarnos solo con los primeros 15
datos_plot <- david_data %>%
    arrange("P-Value") %>% # Ordenar de menor a mayor P-value
    head(15) %>% # Tomar el top 15
    mutate(Term = gsub("^GO:[0-9]+~", "", Term)) # Limpiar el nombre del término GO 

# Genera el Bubble Plot
grafico_burbujas <- ggplot(
    datos_plot,
    aes(x = `Fold Enrichment`, y = reorder(Term, `Fold Enrichment`))
    ) +
    # Crear las burbujas, tamaño por conteo de genes, color por P-Value
    geom_point(aes(size = Count, color = `P-Value`), alpha = 0.8) +
    scale_color_gradient(low = "red", # Rojo para P-values bajos
                        high = "blue") + # Azul para altos
    theme_bw() +
    labs(
        title = "Enriquecimiento de Gene Ontology (Upregulated - CC)",
        x = "Fold Enrichment",
        y = NULL,
        size = "Número de Genes",
        color = "P-value"
    ) +
    theme(
        axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )


# Guarda el gráfico
ggsave(
    "DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/plots/GO_Bubble_Plot_upregulated_CC.png",
    plot = grafico_burbujas,
    width = 11.25,
    height = 7.5,
    dpi = 600
)
