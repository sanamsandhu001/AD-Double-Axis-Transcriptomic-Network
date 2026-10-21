# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 09_GWAS_AD_overlap_robust.R
# Author: Sanampreet Kaur
# Description: Hub definition, GWAS overlap, and Figure 6 generation
# =====================================================

library(dplyr)
library(readr)
library(stringr)
library(igraph)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------
# 1. Load GWAS Catalog File
# -----------------------------------------------------
gwas_path <- "GWAS_AD_associations.tsv"
if (!file.exists(gwas_path)) stop("GWAS_AD_associations.tsv not found in project root.")

cat("\nLoading GWAS catalog...\n")
gwas_raw <- read_tsv(gwas_path, show_col_types = FALSE)

# -----------------------------------------------------
# 2. Filter Alzheimer’s Disease Associations
# -----------------------------------------------------
gwas_sig <- gwas_raw %>%
  dplyr::filter(str_detect(`DISEASE/TRAIT`, regex("Alzheimer", ignore_case = TRUE))) %>%
  dplyr::filter(`P-VALUE` < 5e-8) %>%
  dplyr::distinct(SNP_ID_CURRENT, .keep_all = TRUE)

ad_gene_vector <- gwas_sig$`REPORTED GENE(S)` %>%
  str_split(",|;") %>%
  unlist() %>%
  str_trim() %>%
  toupper() %>%
  unique()
ad_gene_vector <- ad_gene_vector[ad_gene_vector != ""]

cat("Genome-wide significant AD GWAS genes:", length(ad_gene_vector), "\n")

# -----------------------------------------------------
# 3. Load PPI Interaction Networks & Build Graphs
# -----------------------------------------------------
down_edges <- read_csv("data/PPI_Downregulated_interactions.csv", show_col_types = FALSE)
up_edges   <- read_csv("data/PPI_Upregulated_interactions.csv", show_col_types = FALSE)

g_down <- graph_from_data_frame(down_edges, directed = FALSE)
g_up   <- graph_from_data_frame(up_edges, directed = FALSE)

# -----------------------------------------------------
# 4. Compute Centrality Metrics (for Figure 6 plotting only)
# -----------------------------------------------------
# Degree/betweenness for every node are still needed here to
# place all points on the Figure 6 scatter plot. The hub
# classification itself, however, is no longer computed in
# this script — see step 5 below.
down_df <- data.frame(
  Gene = toupper(names(degree(g_down))),
  Degree = as.numeric(degree(g_down)),
  Betweenness = as.numeric(betweenness(g_down, normalized = TRUE)),
  Network = "Downregulated",
  stringsAsFactors = FALSE
)

up_df <- data.frame(
  Gene = toupper(names(degree(g_up))),
  Degree = as.numeric(degree(g_up)),
  Betweenness = as.numeric(betweenness(g_up, normalized = TRUE)),
  Network = "Upregulated",
  stringsAsFactors = FALSE
)

# -----------------------------------------------------
# 5. Load Robust Hubs (defined once, in 05_PPI_network_and_hub_genes.R)
# -----------------------------------------------------
# IMPORTANT: hub identification used to be computed independently
# in this script (top 10% degree + top 10% betweenness, no
# robustness testing) AND in 05_PPI_network_and_hub_genes.R
# (bootstrap edge-resampling on degree only) -- both writing to
# the same results/*_Robust_Hubs.csv files, silently overwriting
# whichever ran second. That's now fixed: 05 is the single
# source of truth for hub identification (dual-centrality
# candidate selection + bootstrap robustness validation, see
# that script for details). This script only *reads* that
# output and must be run after 05.

down_hub_path <- "results/Downregulated_Robust_Hubs.csv"
up_hub_path   <- "results/Upregulated_Robust_Hubs.csv"

if (!file.exists(down_hub_path) || !file.exists(up_hub_path)) {
  stop("Robust hub files not found. Run 05_PPI_network_and_hub_genes.R first.")
}

down_hubs_df <- read_csv(down_hub_path, show_col_types = FALSE)
up_hubs_df   <- read_csv(up_hub_path, show_col_types = FALSE)

down_hubs <- toupper(down_hubs_df$Gene)
up_hubs   <- toupper(up_hubs_df$Gene)

down_df$Is_Hub <- down_df$Gene %in% down_hubs
up_df$Is_Hub   <- up_df$Gene %in% up_hubs

centrality_combined <- dplyr::bind_rows(down_df, up_df) %>%
  dplyr::mutate(
    In_GWAS = Gene %in% ad_gene_vector,
    Category = dplyr::case_when(
      Is_Hub & In_GWAS ~ "Hub & GWAS Risk Gene",
      Is_Hub & !In_GWAS ~ "Robust Network Hub",
      !Is_Hub & In_GWAS ~ "GWAS Gene (Non-Hub)",
      TRUE ~ "Background Node"
    )
  )

down_overlap <- intersect(down_hubs, ad_gene_vector)
up_overlap   <- intersect(up_hubs, ad_gene_vector)

cat("\nRobust Hubs -> Down:", length(down_hubs), "| Up:", length(up_hubs), "\n")
cat("GWAS Overlap -> Down:", paste(down_overlap, collapse=", "), "| Up:", length(up_overlap), "\n")

# -----------------------------------------------------
# 6. Save Tables & Statistical Test
# -----------------------------------------------------
# Note: Robust hub gene lists are owned by 05_PPI_network_and_hub_genes.R
# and are NOT rewritten here anymore (previously this step
# silently overwrote 05's bootstrap-validated output with this
# script's non-bootstrapped quantile version).
write_csv(data.frame(Gene = down_overlap), "results/Downregulated_Hub_GWAS_Overlap.csv")
write_csv(data.frame(Gene = up_overlap), "results/Upregulated_Hub_GWAS_Overlap.csv")
write_csv(centrality_combined, "results/PPI_All_Nodes_Centrality_GWAS.csv")

# -----------------------------------------------------
# 7. FIGURE 6 — Centrality vs. GWAS Overlap Scatter Plot
# -----------------------------------------------------
cat("\nGenerating Figure 6 (Centrality vs GWAS Overlap Scatter Plot)...\n")

centrality_combined$Category <- factor(
  centrality_combined$Category,
  levels = c("Background Node", "GWAS Gene (Non-Hub)", "Robust Network Hub", "Hub & GWAS Risk Gene")
)

# Subset labeled genes (GWAS overlap + top representative hubs).
# Previously this was a hardcoded gene list (including UQCRC1,
# MYD88, STAT5A) that silently went stale once the hub
# definition changed and those genes were no longer robust
# hubs. Now it's derived dynamically from the actual robust hub
# set each time this script runs, so the labels always match
# whatever results/*_Robust_Hubs.csv currently contains: the
# GWAS-overlapping hub(s) plus the top 6 robust hubs by degree
# within each network.
top_hubs_per_network <- centrality_combined %>%
  dplyr::filter(Is_Hub) %>%
  dplyr::group_by(Network) %>%
  dplyr::slice_max(order_by = Degree, n = 6, with_ties = FALSE) %>%
  dplyr::ungroup()

labeled_genes <- centrality_combined %>%
  dplyr::filter(Category == "Hub & GWAS Risk Gene") %>%
  dplyr::bind_rows(top_hubs_per_network) %>%
  dplyr::distinct(Gene, Network, .keep_all = TRUE)

fig6_plot <- ggplot(centrality_combined, aes(x = Degree, y = Betweenness, color = Category)) +
  geom_point(aes(size = Category, alpha = Category)) +
  scale_color_manual(values = c(
    "Background Node" = "grey80",
    "GWAS Gene (Non-Hub)" = "#56B4E9",
    "Robust Network Hub" = "#E69F00",
    "Hub & GWAS Risk Gene" = "#D55E00"
  )) +
  scale_size_manual(values = c(
    "Background Node" = 1.2,
    "GWAS Gene (Non-Hub)" = 2.0,
    "Robust Network Hub" = 2.8,
    "Hub & GWAS Risk Gene" = 4.5
  )) +
  scale_alpha_manual(values = c(
    "Background Node" = 0.35,
    "GWAS Gene (Non-Hub)" = 0.7,
    "Robust Network Hub" = 0.85,
    "Hub & GWAS Risk Gene" = 1.0
  )) +
  geom_text(data = labeled_genes, aes(label = Gene),
            vjust = -0.7, hjust = 0.5, size = 3.5, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Network, scales = "free") +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "PPI Network Centrality and AD GWAS Risk Overlap",
    x = "Degree Centrality (Connectivity)",
    y = "Betweenness Centrality (Normalized)"
  )

ggsave("figures/Figure6_Hub_Centrality_GWAS_Scatter.pdf", fig6_plot, width = 10, height = 5.5)
ggsave("figures/Figure6_Hub_Centrality_GWAS_Scatter.png", fig6_plot, width = 10, height = 5.5, dpi = 300)

cat("Figure 6 saved to figures/Figure6_Hub_Centrality_GWAS_Scatter.png and .pdf\n")
cat("============================================\n")