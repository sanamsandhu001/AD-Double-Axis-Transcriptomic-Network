

# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 03_visualisation_plots.R
# Author: Sanampreet Kaur
# Description: Publication-quality volcano and heatmap
#              (binary-safe RDS loading)
# =====================================================

library(ggplot2)
library(dplyr)
library(pheatmap)

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# Load limma results
# -----------------------------------------------------

limma_results <- read.csv("results/limma_full_results_AD_vs_control.csv")

# -----------------------------------------------------
# Define thresholds
# -----------------------------------------------------

logFC_threshold <- 1
FDR_threshold <- 0.05

# -----------------------------------------------------
# Classify genes
# -----------------------------------------------------

limma_results <- limma_results %>%
  mutate(
    Significance = case_when(
      adj.P.Val < FDR_threshold & logFC > logFC_threshold ~ "Upregulated",
      adj.P.Val < FDR_threshold & logFC < -logFC_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# -----------------------------------------------------
# Volcano Plot
# -----------------------------------------------------

volcano_plot <- ggplot(limma_results,
                       aes(x = logFC,
                           y = -log10(adj.P.Val),
                           color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "firebrick",
    "Downregulated" = "royalblue",
    "Not Significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_threshold),
             linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: AD vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

ggsave("results/plots/Volcano_AD_vs_Control.png",
       volcano_plot,
       width = 8,
       height = 6,
       dpi = 300)

cat("Volcano plot saved.\n")

# -----------------------------------------------------
# Heatmap of Top 50 DEGs
# -----------------------------------------------------

top_genes <- limma_results %>%
  arrange(adj.P.Val) %>%
  slice(1:50) %>%
  pull(GeneSymbol)

# -----------------------------------------------------
# Load expression matrix safely (RDS, not CSV)
# -----------------------------------------------------

expr_matrix <- readRDS("data/expr_GSE5281_clean.rds")

# Ensure rownames exist
if (is.null(rownames(expr_matrix))) {
  stop("Expression matrix does not contain rownames.")
}

# Keep only genes present in expression matrix
top_genes <- intersect(top_genes, rownames(expr_matrix))

heatmap_matrix <- expr_matrix[top_genes, ]

# Safety check
if (nrow(heatmap_matrix) == 0) {
  stop("No matching genes found for heatmap.")
}

pheatmap(
  heatmap_matrix,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  filename = "results/plots/Heatmap_Top50_DEGs.png",
  width = 8,
  height = 10
)

cat("Heatmap saved.\n")

# =====================================================
# END OF SCRIPT 03
# Clean visualization. No binary output.
# =====================================================
