


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 06_Louvain_module_detection.R
# Author: Sanampreet Kaur
# Description: Louvain community detection with
#              stability and random network comparison
# =====================================================

library(igraph)
library(dplyr)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# -----------------------------------------------------
# FUNCTION: Louvain analysis with robustness checks
# -----------------------------------------------------

run_louvain_analysis <- function(input_file, output_prefix) {
  
  cat("\n============================================\n")
  cat("Processing:", input_file, "\n")
  cat("============================================\n")
  
  if (!file.exists(input_file)) {
    stop(paste("Interaction file not found:", input_file))
  }
  
  interactions <- read.csv(input_file, stringsAsFactors = FALSE)
  
  # -----------------------------------------------------
  # Detect edge columns
  # -----------------------------------------------------
  
  if (all(c("protein1_symbol", "protein2_symbol") %in% colnames(interactions))) {
    edge_cols <- c("protein1_symbol", "protein2_symbol")
  } else if (all(c("from", "to") %in% colnames(interactions))) {
    edge_cols <- c("from", "to")
  } else {
    stop("Edge columns not found in interaction file.")
  }
  
  # -----------------------------------------------------
  # Build graph
  # -----------------------------------------------------
  
  g <- igraph::graph_from_data_frame(
    interactions[, edge_cols],
    directed = FALSE
  )
  
  g <- igraph::simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
  
  # -----------------------------------------------------
  # Keep largest connected component
  # -----------------------------------------------------
  
  components_info <- igraph::components(g)
  largest_component <- which.max(components_info$csize)
  
  g <- igraph::induced_subgraph(
    g,
    which(components_info$membership == largest_component)
  )
  
  cat("Nodes in largest component:", vcount(g), "\n")
  cat("Edges:", ecount(g), "\n")
  
  # -----------------------------------------------------
  # Louvain clustering (real network)
  # -----------------------------------------------------
  
  set.seed(123)
  louvain_result <- igraph::cluster_louvain(g)
  
  modularity_real <- igraph::modularity(louvain_result)
  cat("Real network modularity:", modularity_real, "\n")
  
  membership_df <- data.frame(
    GeneSymbol = names(igraph::membership(louvain_result)),
    Module = as.numeric(igraph::membership(louvain_result)),
    stringsAsFactors = FALSE
  )
  
  module_summary <- membership_df %>%
    group_by(Module) %>%
    summarise(Size = n(), .groups = "drop") %>%
    arrange(desc(Size))
  
  print(module_summary)
  
  # -----------------------------------------------------
  # Stability check (10 runs)
  # -----------------------------------------------------
  
  modularity_runs <- numeric(10)
  
  for (i in 1:10) {
    set.seed(i * 100)
    temp_cluster <- igraph::cluster_louvain(g)
    modularity_runs[i] <- igraph::modularity(temp_cluster)
  }
  
  mean_mod <- mean(modularity_runs)
  sd_mod   <- sd(modularity_runs)
  
  cat("Mean modularity (10 runs):", mean_mod, "\n")
  cat("SD modularity (10 runs):", sd_mod, "\n")
  
  # -----------------------------------------------------
  # Random network comparison
  # -----------------------------------------------------
  
  cat("\nGenerating random network for comparison...\n")
  
  set.seed(999)
  
  g_random <- igraph::sample_gnm(
    n = vcount(g),
    m = ecount(g),
    directed = FALSE,
    loops = FALSE
  )
  
  random_cluster <- igraph::cluster_louvain(g_random)
  modularity_random <- igraph::modularity(random_cluster)
  
  cat("Random network modularity:", modularity_random, "\n")
  
  # -----------------------------------------------------
  # Save results
  # -----------------------------------------------------
  
  write.csv(
    membership_df,
    file.path(results_dir,
              paste0("Louvain_", output_prefix,
                     "_Module_Assignments.csv")),
    row.names = FALSE
  )
  
  write.csv(
    module_summary,
    file.path(results_dir,
              paste0("Louvain_", output_prefix,
                     "_Module_Summary.csv")),
    row.names = FALSE
  )
  
  robustness_stats <- data.frame(
    Real_Modularity = modularity_real,
    Mean_Modularity_10_Runs = mean_mod,
    SD_Modularity_10_Runs = sd_mod,
    Random_Modularity = modularity_random
  )
  
  write.csv(
    robustness_stats,
    file.path(results_dir,
              paste0("Louvain_", output_prefix,
                     "_Robustness_Stats.csv")),
    row.names = FALSE
  )
  
  cat("\nRobustness statistics saved.\n")
  cat("============================================\n")
}

# -----------------------------------------------------
# RUN ANALYSIS
# -----------------------------------------------------

run_louvain_analysis(
  "data/PPI_Downregulated_interactions.csv",
  "Downregulated"
)

run_louvain_analysis(
  "data/PPI_Upregulated_interactions.csv",
  "Upregulated"
)

# =====================================================
# END OF SCRIPT 06
# Louvain clustering with stability and null comparison
# =====================================================
  