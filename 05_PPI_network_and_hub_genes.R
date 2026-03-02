


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 05_PPI_network_and_hub_genes.R
# Author: Sanampreet Kaur
# Description: STRING PPI network construction,
#              filtering, and topology metrics
# =====================================================

library(STRINGdb)
library(dplyr)
library(igraph)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

options(timeout = 300)

# -----------------------------------------------------
# Initialize STRING
# -----------------------------------------------------

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 700,
  input_directory = ""
)

# -----------------------------------------------------
# FUNCTION: Build PPI network
# -----------------------------------------------------

build_ppi_network <- function(deg_file, output_prefix) {
  
  cat("\n============================================\n")
  cat("Processing:", deg_file, "\n")
  cat("============================================\n")
  
  if (!file.exists(deg_file)) {
    stop(paste("DEG file not found:", deg_file))
  }
  
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)
  deg$GeneSymbol <- toupper(trimws(deg$GeneSymbol))
  
  # ---------------- Map to STRING ----------------
  
  deg_mapped <- string_db$map(
    deg,
    "GeneSymbol",
    removeUnmappedRows = TRUE
  )
  
  cat("Input genes:", nrow(deg), "\n")
  cat("Mapped to STRING:", nrow(deg_mapped), "\n")
  
  if (nrow(deg_mapped) < 10) {
    stop("Too few genes mapped to STRING.")
  }
  
  # ---------------- Get interactions ----------------
  
  interactions_raw <- string_db$get_interactions(deg_mapped$STRING_id)
  
  interactions_filtered <- interactions_raw %>%
    filter(combined_score >= 700) %>%
    filter(from %in% deg_mapped$STRING_id &
             to   %in% deg_mapped$STRING_id) %>%
    distinct(from, to)
  
  # ---------------- Convert STRING IDs → Gene Symbols ----------------
  
  id_map <- deg_mapped %>%
    select(STRING_id, GeneSymbol) %>%
    distinct()
  
  interactions_mapped <- interactions_filtered %>%
    left_join(id_map, by = c("from" = "STRING_id")) %>%
    rename(protein1_symbol = GeneSymbol) %>%
    left_join(id_map, by = c("to" = "STRING_id")) %>%
    rename(protein2_symbol = GeneSymbol) %>%
    filter(!is.na(protein1_symbol) &
             !is.na(protein2_symbol)) %>%
    distinct(protein1_symbol, protein2_symbol)
  
  if (nrow(interactions_mapped) == 0) {
    stop("No valid interactions after symbol conversion.")
  }
  
  # ---------------- Build graph ----------------
  
  g <- graph_from_data_frame(
    interactions_mapped[, c("protein1_symbol", "protein2_symbol")],
    directed = FALSE
  )
  
  g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
  
  # ---------------- Largest connected component ----------------
  
  components_info <- components(g)
  largest_component <- which.max(components_info$csize)
  g <- induced_subgraph(
    g,
    which(components_info$membership == largest_component)
  )
  
  cat("Nodes in largest component:", vcount(g), "\n")
  cat("Edges:", ecount(g), "\n")
  
  # ---------------- Save interactions for Louvain ----------------
  
  edge_matrix <- igraph::ends(g, E(g))
  
  interactions_df <- data.frame(
    from = edge_matrix[, 1],
    to   = edge_matrix[, 2],
    stringsAsFactors = FALSE
  )
  
  write.csv(
    interactions_df,
    paste0("data/PPI_", output_prefix, "_interactions.csv"),
    row.names = FALSE
  )
  
  # ---------------- Network metrics ----------------
  
  network_metrics <- data.frame(
    Nodes = vcount(g),
    Edges = ecount(g),
    Average_Degree = mean(degree(g)),
    Clustering_Coefficient = transitivity(g, type = "average"),
    Density = edge_density(g)
  )
  
  write.csv(
    network_metrics,
    paste0("results/PPI_", output_prefix, "_network_metrics.csv"),
    row.names = FALSE
  )
  
  cat("Network metrics saved.\n")
  cat("============================================\n")
}

# -----------------------------------------------------
# RUN FOR BOTH DIRECTIONS
# -----------------------------------------------------

build_ppi_network(
  "results/DEG_downregulated_AD_vs_control.csv",
  "Downregulated"
)

build_ppi_network(
  "results/DEG_upregulated_AD_vs_control.csv",
  "Upregulated"
)

# =====================================================
# END OF SCRIPT 05
# Clean STRING PPI networks constructed
# =====================================================
