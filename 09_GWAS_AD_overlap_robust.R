


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 09_GWAS_AD_overlap.R
# Author: Sanampreet Kaur
# Description: Robust hub definition and overlap with
#              genome-wide significant AD GWAS genes
# =====================================================

library(dplyr)
library(readr)
library(stringr)
library(igraph)

dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------
# 1. Load GWAS Catalog File
# -----------------------------------------------------

gwas_path <- "GWAS_AD_associations.tsv"

if (!file.exists(gwas_path)) {
  stop("GWAS_AD_associations.tsv not found in project root.")
}

cat("\nLoading GWAS catalog...\n")

gwas_raw <- read_tsv(gwas_path, show_col_types = FALSE)

# -----------------------------------------------------
# 2. Filter Alzheimer’s Disease Associations
#    + Genome-wide significance
#    + Unique lead SNPs
# -----------------------------------------------------

gwas_sig <- gwas_raw %>%
  filter(str_detect(`DISEASE/TRAIT`,
                    regex("Alzheimer", ignore_case = TRUE))) %>%
  filter(`P-VALUE` < 5e-8) %>%
  distinct(SNP_ID_CURRENT, .keep_all = TRUE)

# -----------------------------------------------------
# 3. Extract REPORTED Genes
# -----------------------------------------------------

ad_gene_vector <- gwas_sig$`REPORTED GENE(S)` %>%
  str_split(",|;") %>%
  unlist() %>%
  str_trim() %>%
  toupper() %>%
  unique()

ad_gene_vector <- ad_gene_vector[ad_gene_vector != ""]

cat("Genome-wide significant AD GWAS genes:",
    length(ad_gene_vector), "\n")

# -----------------------------------------------------
# 4. Load PPI Interaction Networks
# -----------------------------------------------------

down_edges_path <- "data/PPI_Downregulated_interactions.csv"
up_edges_path   <- "data/PPI_Upregulated_interactions.csv"

if (!file.exists(down_edges_path) | !file.exists(up_edges_path)) {
  stop("PPI interaction files not found. Run Script 05 first.")
}

down_edges <- read_csv(down_edges_path, show_col_types = FALSE)
up_edges   <- read_csv(up_edges_path, show_col_types = FALSE)

# -----------------------------------------------------
# 5. Build Graph Objects
# -----------------------------------------------------

g_down <- graph_from_data_frame(down_edges, directed = FALSE)
g_up   <- graph_from_data_frame(up_edges, directed = FALSE)

# -----------------------------------------------------
# 6. Compute Centrality Metrics
# -----------------------------------------------------

down_df <- data.frame(
  Gene = toupper(names(degree(g_down))),
  Degree = as.numeric(degree(g_down)),
  Betweenness = as.numeric(betweenness(g_down, normalized = TRUE))
)

up_df <- data.frame(
  Gene = toupper(names(degree(g_up))),
  Degree = as.numeric(degree(g_up)),
  Betweenness = as.numeric(betweenness(g_up, normalized = TRUE))
)

# -----------------------------------------------------
# 7. Define Robust Hubs (Top 10% Degree + Betweenness)
# -----------------------------------------------------

deg_thresh_down <- quantile(down_df$Degree, 0.9)
bet_thresh_down <- quantile(down_df$Betweenness, 0.9)

down_hubs <- down_df %>%
  filter(Degree >= deg_thresh_down &
           Betweenness >= bet_thresh_down) %>%
  pull(Gene)

deg_thresh_up <- quantile(up_df$Degree, 0.9)
bet_thresh_up <- quantile(up_df$Betweenness, 0.9)

up_hubs <- up_df %>%
  filter(Degree >= deg_thresh_up &
           Betweenness >= bet_thresh_up) %>%
  pull(Gene)

cat("\nRobust hubs identified:\n")
cat("Downregulated:", length(down_hubs), "\n")
cat("Upregulated:", length(up_hubs), "\n")

# -----------------------------------------------------
# 8. Overlap with GWAS Genes
# -----------------------------------------------------

down_overlap <- intersect(down_hubs, ad_gene_vector)
up_overlap   <- intersect(up_hubs, ad_gene_vector)

cat("\nHub–GWAS overlap:\n")
cat("Downregulated overlap:", length(down_overlap), "\n")
cat("Upregulated overlap:", length(up_overlap), "\n")

# -----------------------------------------------------
# 9. Fisher’s Exact Test for Enrichment
# -----------------------------------------------------

all_down_genes <- down_df$Gene
all_up_genes   <- up_df$Gene

matrix_down <- matrix(c(
  length(down_overlap),
  length(down_hubs) - length(down_overlap),
  sum(all_down_genes %in% ad_gene_vector) - length(down_overlap),
  length(all_down_genes) -
    sum(all_down_genes %in% ad_gene_vector) -
    (length(down_hubs) - length(down_overlap))
), nrow = 2)

matrix_up <- matrix(c(
  length(up_overlap),
  length(up_hubs) - length(up_overlap),
  sum(all_up_genes %in% ad_gene_vector) - length(up_overlap),
  length(all_up_genes) -
    sum(all_up_genes %in% ad_gene_vector) -
    (length(up_hubs) - length(up_overlap))
), nrow = 2)

fisher_down <- fisher.test(matrix_down)
fisher_up   <- fisher.test(matrix_up)

cat("\nFisher Test Results:\n")
cat("Downregulated p-value:", fisher_down$p.value, "\n")
cat("Upregulated p-value:", fisher_up$p.value, "\n")

# -----------------------------------------------------
# 10. Save Outputs
# -----------------------------------------------------

write_csv(data.frame(Gene = down_hubs),
          "results/Downregulated_Robust_Hubs.csv")

write_csv(data.frame(Gene = up_hubs),
          "results/Upregulated_Robust_Hubs.csv")

write_csv(data.frame(Gene = down_overlap),
          "results/Downregulated_Hub_GWAS_Overlap.csv")

write_csv(data.frame(Gene = up_overlap),
          "results/Upregulated_Hub_GWAS_Overlap.csv")

cat("\nGWAS overlap analysis complete.\n")

# =====================================================
# END OF SCRIPT 09
# Robust GWAS–Network integration complete
# =====================================================