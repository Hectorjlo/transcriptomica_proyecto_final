# Load libraries
library(ggplot2)
library(stringr)

# Read the results file
path <- "results/stats/final_stats.txt"
data_text <- paste(readLines(path), collapse = "\n")

# Parse text lines
lines <- readLines(path)

# Initialize all the variables
assemblers <- c()
srrs <- c()
types <- c()
overall <- c()
non_aligned <- c()
once <- c()
multi <- c()

for (i in seq_along(lines)) {
  # Parse the "fasta-like" file to the previous variables
  if (grepl("^>ASSEMBLER:", lines[i])) {
    # Assign them trough REGEX patterns
    assemblers <- c(assemblers, str_extract(lines[i], "(?<=ASSEMBLER: )\\w+"))
    srrs <- c(srrs, str_extract(lines[i], "SRR\\d+"))
    types <- c(types, str_extract(lines[i], "(?<=\\()PE|SE(?=\\))"))
    overall <- c(overall, as.numeric(str_extract(lines[i + 1], "[\\d.]+(?=%)")))
    non_aligned <- c(
      non_aligned,
      as.numeric(str_extract(lines[i + 2], "[\\d.]+(?=%)"))
    )
    once <- c(once, as.numeric(str_extract(lines[i + 3], "[\\d.]+(?=%)")))
    multi <- c(multi, as.numeric(str_extract(lines[i + 4], "[\\d.]+")))
  }
}

# Create a data.frame
df_wide <- data.frame(
  Assembler = assemblers,
  SRR = srrs,
  Type = types,
  Overall = overall,
  Non_Aligned = non_aligned,
  Once = once,
  Multi = multi
)

# Transform the data.frame to a long format, easier to manipulate
# It expands rows to columns to build a longer data.frame
df_long <- tidyr::pivot_longer(
  df_wide,
  cols = c(Overall, Non_Aligned, Once, Multi),
  names_to = "Metric",
  values_to = "Percentage"
)

# Create a new column, with the Assembler and Type in the same field
df_long$Group <- paste(df_long$Assembler, df_long$Type, sep = " - ")

# Tags for the metrics
df_long$Metric <- factor(
  df_long$Metric,
  levels = c("Overall", "Non_Aligned", "Once", "Multi"),
  labels = c(
    "Overall Alignment",
    "Non Aligned",
    "Aligned Once",
    "Aligned >1 times"
  )
)

# Geom point plot, showing of stat with differece in color and shape
point_plot <- ggplot(
  df_long,
  aes(x = SRR, y = Percentage, color = Group, shape = Metric)
) +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Alignment Statistics per SRR",
    subtitle = "All metrics by Assembler and Read Type (PE vs SE)",
    x = "SRR Accession",
    y = "Percentage (%)",
    color = "Condition",
    shape = "Metric"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 100))

# Save in high quality the plot
ggsave(
  "results/plots/stats_points.png",
  plot = point_plot,
  dpi = 800,
  width = 11.25,
  height = 7.5,
  bg = "white"
)
