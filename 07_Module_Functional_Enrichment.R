


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 07_Module_Functional_Enrichment.R
# Author: Sanampreet Kaur
# Description: Functional enrichment of Louvain modules
#              with dataset-specific background
# =====================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

results_dir <- file.path(getwd(), "results/module_enrichment")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# DEFINE BACKGROUND UNIVERSE
# -----------------------------------------------------

cat("\nDefining background universe from limma results...\n")

limma_full <- read.csv("data/limma_full_results_AD_vs_control.csv")

limma_full$GeneSymbol <- toupper(trimws(limma_full$GeneSymbol))

all_symbols <- unique(limma_full$GeneSymbol)

all_entrez <- clusterProfiler::bitr(
  all_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

background_universe <- unique(all_entrez$ENTREZID)

cat("Background Entrez IDs:", length(background_universe), "\n")

# -----------------------------------------------------
# FUNCTION: Enrich Louvain modules
# -----------------------------------------------------

run_module_enrichment <- function(module_file, network_label, top_n = 5) {
  
  cat("\n============================================\n")
  cat("Processing:", network_label, "\n")
  cat("============================================\n")
  
  modules <- read.csv(module_file)
  
  modules$GeneSymbol <- toupper(trimws(modules$GeneSymbol))
  
  # Get largest modules
  module_sizes <- modules %>%
    group_by(Module) %>%
    summarise(Size = n(), .groups = "drop") %>%
    arrange(desc(Size))
  
  print(module_sizes)
  
  top_modules <- module_sizes$Module[1:min(top_n, nrow(module_sizes))]
  
  cat("Top modules selected:", top_modules, "\n")
  
  for (mod in top_modules) {
    
    cat("\n--- Module", mod, "---\n")
    
    gene_symbols <- modules %>%
      filter(Module == mod) %>%
      pull(GeneSymbol)
    
    if (length(gene_symbols) < 15) {
      cat("Module too small. Skipping.\n")
      next
    }
    
    gene_entrez <- clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
    
    gene_list <- unique(gene_entrez$ENTREZID)
    
    cat("Genes in module:", length(gene_symbols), "\n")
    cat("Mapped Entrez IDs:", length(gene_list), "\n")
    
    if (length(gene_list) < 10) {
      cat("Too few mapped genes. Skipping.\n")
      next
    }
    
    # ---------------- GO BP ----------------
    
    ego <- clusterProfiler::enrichGO(
      gene          = gene_list,
      universe      = background_universe,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      
      ego_simplified <- clusterProfiler::simplify(
        ego,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min
      )
      
      write.csv(
        as.data.frame(ego_simplified),
        file.path(results_dir,
                  paste0(network_label,
                         "_Module_", mod,
                         "_GO_BP.csv")),
        row.names = FALSE
      )
      
      cat("GO terms:", nrow(as.data.frame(ego_simplified)), "\n")
      
    } else {
      cat("No significant GO terms.\n")
    }
    
    # ---------------- KEGG ----------------
    
    ekegg <- clusterProfiler::enrichKEGG(
      gene         = gene_list,
      universe     = background_universe,
      organism     = "hsa",
      pvalueCutoff = 0.05
    )
    
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      
      write.csv(
        as.data.frame(ekegg),
        file.path(results_dir,
                  paste0(network_label,
                         "_Module_", mod,
                         "_KEGG.csv")),
        row.names = FALSE
      )
      
      cat("KEGG pathways:", nrow(as.data.frame(ekegg)), "\n")
      
    } else {
      cat("No significant KEGG pathways.\n")
    }
  }
  
  cat("\nModule enrichment completed for", network_label, "\n")
  cat("============================================\n")
}

# -----------------------------------------------------
# RUN FOR BOTH NETWORKS
# -----------------------------------------------------

run_module_enrichment(
  "results/Louvain_Downregulated_Module_Assignments.csv",
  "Downregulated",
  top_n = 5
)

run_module_enrichment(
  "results/Louvain_Upregulated_Module_Assignments.csv",
  "Upregulated",
  top_n = 5
)

# =====================================================
# END OF SCRIPT 07
# Network-structured module enrichment complete
# =====================================================




# =====================================================
# Project: Molecular network characterization of Alzheimer’s disease
# Script: 07_Module_Functional_Enrichment.R
# Description: Functional enrichment of Louvain modules
#              with dataset-specific background
# =====================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

results_dir <- "results/module_enrichment"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# DEFINE BACKGROUND UNIVERSE
# -----------------------------------------------------

cat("\nDefining background universe from limma results...\n")

limma_path <- "results/limma_full_results_AD_vs_control.csv"

if (!file.exists(limma_path)) {
  stop("Limma results not found. Run Script 02 first.")
}

limma_full <- read.csv(limma_path, stringsAsFactors = FALSE)
limma_full$GeneSymbol <- toupper(trimws(limma_full$GeneSymbol))

all_symbols <- unique(limma_full$GeneSymbol)

all_entrez <- clusterProfiler::bitr(
  all_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

background_universe <- unique(all_entrez$ENTREZID)

cat("Background Entrez IDs:", length(background_universe), "\n")

# -----------------------------------------------------
# FUNCTION: Enrich Louvain modules
# -----------------------------------------------------

run_module_enrichment <- function(module_file, network_label, top_n = 5) {
  
  cat("\n============================================\n")
  cat("Processing:", network_label, "\n")
  cat("============================================\n")
  
  if (!file.exists(module_file)) {
    stop(paste("Module file not found:", module_file))
  }
  
  modules <- read.csv(module_file, stringsAsFactors = FALSE)
  modules$GeneSymbol <- toupper(trimws(modules$GeneSymbol))
  
  module_sizes <- modules %>%
    group_by(Module) %>%
    summarise(Size = n(), .groups = "drop") %>%
    arrange(desc(Size))
  
  print(module_sizes)
  
  top_modules <- module_sizes$Module[1:min(top_n, nrow(module_sizes))]
  cat("Top modules selected:", top_modules, "\n")
  
  for (mod in top_modules) {
    
    cat("\n--- Module", mod, "---\n")
    
    gene_symbols <- modules %>%
      filter(Module == mod) %>%
      pull(GeneSymbol)
    
    if (length(gene_symbols) < 15) {
      cat("Module too small. Skipping.\n")
      next
    }
    
    gene_entrez <- clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
    
    gene_list <- unique(gene_entrez$ENTREZID)
    
    cat("Genes in module:", length(gene_symbols), "\n")
    cat("Mapped Entrez IDs:", length(gene_list), "\n")
    
    if (length(gene_list) < 10) {
      cat("Too few mapped genes. Skipping.\n")
      next
    }
    
    # ---------------- GO BP ----------------
    
    ego <- clusterProfiler::enrichGO(
      gene          = gene_list,
      universe      = background_universe,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      
      ego_simplified <- clusterProfiler::simplify(
        ego,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min
      )
      
      write.csv(
        as.data.frame(ego_simplified),
        file.path(results_dir,
                  paste0(network_label,
                         "_Module_", mod,
                         "_GO_BP.csv")),
        row.names = FALSE
      )
      
      cat("GO terms:", nrow(as.data.frame(ego_simplified)), "\n")
      
    } else {
      cat("No significant GO terms.\n")
    }
    
    # ---------------- KEGG ----------------
    
    ekegg <- clusterProfiler::enrichKEGG(
      gene         = gene_list,
      universe     = background_universe,
      organism     = "hsa",
      pvalueCutoff = 0.05
    )
    
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      
      write.csv(
        as.data.frame(ekegg),
        file.path(results_dir,
                  paste0(network_label,
                         "_Module_", mod,
                         "_KEGG.csv")),
        row.names = FALSE
      )
      
      cat("KEGG pathways:", nrow(as.data.frame(ekegg)), "\n")
      
    } else {
      cat("No significant KEGG pathways.\n")
    }
  }
  
  cat("\nModule enrichment completed for", network_label, "\n")
  cat("============================================\n")
}

# -----------------------------------------------------
# RUN FOR BOTH NETWORKS
# -----------------------------------------------------

run_module_enrichment(
  "results/Louvain_Downregulated_Module_Assignments.csv",
  "Downregulated",
  top_n = 5
)

run_module_enrichment(
  "results/Louvain_Upregulated_Module_Assignments.csv",
  "Upregulated",
  top_n = 5
)

# =====================================================
# END OF SCRIPT 07
# Network-structured module enrichment complete
# =====================================================