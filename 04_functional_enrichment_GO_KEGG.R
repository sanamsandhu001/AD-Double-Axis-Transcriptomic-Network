


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 04_functional_enrichment_GO_KEGG.R
# Author: Sanampreet Kaur
# Description: Whole-DEG enrichment (FDR-based)
#              with dataset-specific background
# =====================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# DEFINE BACKGROUND UNIVERSE
# -----------------------------------------------------

cat("\nDefining background universe from limma results...\n")

limma_path <- "results/limma_full_results_AD_vs_control.csv"

if (!file.exists(limma_path)) {
  stop("Limma full results not found in results/. Run Script 02 first.")
}

limma_full <- read.csv(limma_path)

# Clean gene symbols
limma_full$GeneSymbol <- toupper(trimws(limma_full$GeneSymbol))

all_symbols <- unique(limma_full$GeneSymbol)

all_entrez <- clusterProfiler::bitr(
  all_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

background_universe <- unique(all_entrez$ENTREZID)

cat("Total genes tested:", length(all_symbols), "\n")
cat("Mapped background Entrez IDs:", length(background_universe), "\n")

# -----------------------------------------------------
# FUNCTION: Run enrichment safely
# -----------------------------------------------------

run_enrichment <- function(deg_file, label) {
  
  cat("\n============================================\n")
  cat("Processing:", label, "\n")
  cat("============================================\n")
  
  if (!file.exists(deg_file)) {
    stop(paste("DEG file not found:", deg_file))
  }
  
  deg <- read.csv(deg_file)
  
  if (nrow(deg) < 10) {
    cat("Too few genes for enrichment. Skipping.\n")
    return(NULL)
  }
  
  # Clean symbols
  deg$GeneSymbol <- toupper(trimws(deg$GeneSymbol))
  
  gene_entrez <- clusterProfiler::bitr(
    deg$GeneSymbol,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  gene_list <- unique(gene_entrez$ENTREZID)
  
  cat("Input genes:", nrow(deg), "\n")
  cat("Mapped Entrez IDs:", length(gene_list), "\n")
  
  if (length(gene_list) < 10) {
    cat("Too few mapped genes for enrichment. Skipping.\n")
    return(NULL)
  }
  
  # ---------------- GO Biological Process ----------------
  
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
      paste0("results/enrichment/",
             label, "_GO_BP_universe_corrected.csv"),
      row.names = FALSE
    )
    
    cat("GO terms saved:", nrow(as.data.frame(ego_simplified)), "\n")
  } else {
    cat("No significant GO terms found.\n")
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
      paste0("results/enrichment/",
             label, "_KEGG_universe_corrected.csv"),
      row.names = FALSE
    )
    
    cat("KEGG pathways saved:", nrow(as.data.frame(ekegg)), "\n")
  } else {
    cat("No significant KEGG pathways found.\n")
  }
  
  cat("Enrichment completed for", label, "\n")
  cat("============================================\n")
}

# -----------------------------------------------------
# RUN FOR BOTH DIRECTIONS
# -----------------------------------------------------

run_enrichment(
  "results/DEG_downregulated_AD_vs_control.csv",
  "Downregulated"
)

run_enrichment(
  "results/DEG_upregulated_AD_vs_control.csv",
  "Upregulated"
)

# =====================================================
# END OF SCRIPT 04
# Dataset-specific, redundancy-reduced enrichment complete
# =====================================================