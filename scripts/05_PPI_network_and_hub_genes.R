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
    dplyr::filter(combined_score >= 700) %>%
    dplyr::filter(from %in% deg_mapped$STRING_id &
                    to   %in% deg_mapped$STRING_id) %>%
    dplyr::distinct(from, to)
  
  # ---------------- Convert STRING IDs → Gene Symbols ----------------
  
  id_map <- deg_mapped %>%
    dplyr::select(STRING_id, GeneSymbol) %>%
    dplyr::distinct()
  
  interactions_mapped <- interactions_filtered %>%
    dplyr::left_join(id_map, by = c("from" = "STRING_id")) %>%
    dplyr::rename(protein1_symbol = GeneSymbol) %>%
    dplyr::left_join(id_map, by = c("to" = "STRING_id")) %>%
    dplyr::rename(protein2_symbol = GeneSymbol) %>%
    dplyr::filter(!is.na(protein1_symbol) &
                    !is.na(protein2_symbol)) %>%
    dplyr::distinct(protein1_symbol, protein2_symbol)
  
  if (nrow(interactions_mapped) == 0) {
    stop("No valid interactions after symbol conversion.")
  }
  
  # ---------------- Build graph ----------------
  
  g <- graph_from_data_frame(
    interactions_mapped[, c("protein1_symbol", "protein2_symbol")],
    directed = FALSE
  )
  
  g <- igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
  
  # ---------------- Largest connected component ----------------
  
  components_info <- igraph::components(g)
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
# =====================================================
# ADDITION TO: 05_PPI_network_and_hub_genes.R
# Section: Robust hub gene identification
# -----------------------------------------------------
# This is the single, canonical hub-gene definition for the
# project. It replaces two earlier, inconsistent versions that
# both wrote to results/*_Robust_Hubs.csv (a degree-only
# bootstrap version that lived here, and a dual-centrality,
# non-bootstrapped version that lived in
# 09_GWAS_AD_overlap_robust.R). Script 09 now reads its hub
# list from the output of this function instead of recomputing
# its own — see the updated 09_GWAS_AD_overlap_robust.R.
#
# Method (two stages):
#   Stage 1 - Candidate selection: compute degree centrality and
#     normalized betweenness centrality on the full network. A
#     gene is a *candidate* hub if it falls in the top decile
#     (>= 90th percentile) of BOTH measures simultaneously.
#     This is proportional (scales with network size) rather
#     than a fixed top-N cutoff, so the downregulated and
#     upregulated networks — which differ >2x in node count —
#     are held to the same relative standard.
#   Stage 2 - Bootstrap robustness validation: repeat
#     n_iterations times: randomly remove a fraction
#     (edge_dropout) of edges (all nodes retained), recompute
#     degree, betweenness, and their 90th-percentile thresholds
#     on the perturbed network, and record whether each
#     candidate still meets the dual top-decile criterion in
#     that iteration. Retain as a final "robust hub gene" only
#     candidates that hit in at least robustness_threshold
#     proportion of iterations.
#
# These parameters (percentile cutoff, edge_dropout,
# n_iterations, robustness_threshold) are analytical choices —
# sanity-check them for your own data before treating the
# results as final.
# =====================================================

library(igraph)
library(dplyr)

identify_hub_genes <- function(interactions_file,
                               output_prefix,
                               percentile = 0.90,
                               edge_dropout = 0.10,
                               n_iterations = 100,
                               robustness_threshold = 0.90,
                               seed = 42) {
  
  cat("\n============================================\n")
  cat("Identifying hub genes:", output_prefix, "\n")
  cat("============================================\n")
  
  set.seed(seed)
  
  interactions_df <- read.csv(interactions_file, stringsAsFactors = FALSE)
  
  g <- graph_from_data_frame(
    interactions_df[, c("from", "to")],
    directed = FALSE
  )
  g <- igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
  
  cat("Network loaded: ", vcount(g), "nodes,", ecount(g), "edges\n")
  
  # ---------------- Stage 1: dual-centrality candidate hubs ----------------
  
  full_degree      <- degree(g)
  full_betweenness <- betweenness(g, normalized = TRUE)
  
  deg_thresh <- quantile(full_degree, percentile)
  bet_thresh <- quantile(full_betweenness, percentile)
  
  is_candidate <- (full_degree >= deg_thresh) & (full_betweenness >= bet_thresh)
  candidate_hubs <- names(full_degree)[is_candidate]
  
  cat("Candidate hubs (top", (1 - percentile) * 100,
      "% degree AND betweenness):", length(candidate_hubs), "\n")
  cat("  degree threshold:", round(deg_thresh, 3),
      " betweenness threshold:", round(bet_thresh, 5), "\n")
  
  # ---------------- Stage 2: bootstrap robustness filtering ----------------
  
  hit_count <- setNames(rep(0, length(candidate_hubs)), candidate_hubs)
  n_edges <- ecount(g)
  keep_n  <- round(n_edges * (1 - edge_dropout))
  
  for (i in seq_len(n_iterations)) {
    
    keep_idx <- sample(seq_len(n_edges), size = keep_n)
    g_perturbed <- subgraph.edges(g, keep_idx, delete.vertices = FALSE)
    
    perturbed_degree      <- degree(g_perturbed)
    perturbed_betweenness <- betweenness(g_perturbed, normalized = TRUE)
    
    deg_thresh_iter <- quantile(perturbed_degree, percentile)
    bet_thresh_iter <- quantile(perturbed_betweenness, percentile)
    
    is_hub_iter <- (perturbed_degree[candidate_hubs] >= deg_thresh_iter) &
      (perturbed_betweenness[candidate_hubs] >= bet_thresh_iter)
    
    hits <- candidate_hubs[is_hub_iter]
    hit_count[hits] <- hit_count[hits] + 1
  }
  
  robustness_score <- hit_count / n_iterations
  
  robust_hubs <- names(robustness_score)[robustness_score >= robustness_threshold]
  
  cat("Robust hubs (>=", robustness_threshold * 100, "% of iterations):",
      length(robust_hubs), "out of", length(candidate_hubs), "candidates\n")
  
  # ---------------- Save results ----------------
  # Full candidate table (useful for auditing/adjusting thresholds later)
  
  candidate_df <- data.frame(
    Gene = candidate_hubs,
    Degree = full_degree[candidate_hubs],
    Betweenness = full_betweenness[candidate_hubs],
    RobustnessScore = robustness_score[candidate_hubs],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(desc(RobustnessScore), desc(Degree))
  
  write.csv(
    candidate_df,
    paste0("results/", output_prefix, "_Hub_Candidates_All.csv"),
    row.names = FALSE
  )
  
  # Final robust hub set only
  
  hub_df <- candidate_df %>%
    dplyr::filter(RobustnessScore >= robustness_threshold)
  
  write.csv(
    hub_df,
    paste0("results/", output_prefix, "_Robust_Hubs.csv"),
    row.names = FALSE
  )
  
  cat("Saved to results/", output_prefix, "_Robust_Hubs.csv\n", sep = "")
  cat("Saved to results/", output_prefix, "_Hub_Candidates_All.csv\n", sep = "")
  cat("============================================\n")
  
  return(hub_df)
}

# -----------------------------------------------------
# RUN FOR BOTH DIRECTIONS
# -----------------------------------------------------
# Reads the edge lists already written by build_ppi_network()
# earlier in this script.

hubs_down <- identify_hub_genes(
  "data/PPI_Downregulated_interactions.csv",
  "Downregulated"
)

hubs_up <- identify_hub_genes(
  "data/PPI_Upregulated_interactions.csv",
  "Upregulated"
)

# =====================================================
# END OF ADDITION
# Robust hub gene identification complete
# =====================================================