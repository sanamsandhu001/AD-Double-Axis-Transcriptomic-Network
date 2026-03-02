


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 02_differential_expression_limma.R
# Author: Sanampreet Kaur
# Description: Differential expression using limma
#              FDR-based directional filtering
# =====================================================

library(limma)
library(dplyr)

dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------
# Load cleaned expression and phenotype data
# -----------------------------------------------------

expr_matrix <- readRDS("data/expr_GSE5281_clean.rds")
pheno_data  <- readRDS("data/pdata_GSE5281_clean.rds")

expr_matrix <- as.matrix(expr_matrix)

# -----------------------------------------------------
# Strict alignment (no intersect)
# -----------------------------------------------------

pheno_data <- pheno_data[colnames(expr_matrix), ]
stopifnot(all(colnames(expr_matrix) == rownames(pheno_data)))

cat("Samples aligned:", ncol(expr_matrix), "\n")

# -----------------------------------------------------
# Define group variable
# -----------------------------------------------------

pheno_data$Group <- factor(pheno_data$Group,
                           levels = c("Control", "AD"))

# -----------------------------------------------------
# Design matrix
# -----------------------------------------------------

design <- model.matrix(~ Group, data = pheno_data)
colnames(design) <- make.names(colnames(design))

# -----------------------------------------------------
# Fit linear model
# -----------------------------------------------------

fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

# -----------------------------------------------------
# Extract results
# -----------------------------------------------------

results <- topTable(
  fit,
  coef = "GroupAD",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# -----------------------------------------------------
# FIX: Attach correct gene symbols
# -----------------------------------------------------

results$GeneSymbol <- rownames(expr_matrix)

# -----------------------------------------------------
# Save full results
# -----------------------------------------------------

write.csv(
  results,
  "results/limma_full_results_AD_vs_control.csv",
  row.names = FALSE
)

# -----------------------------------------------------
# Define FDR significance (NO change to logic)
# -----------------------------------------------------

FDR_threshold <- 0.05

deg_significant <- results %>%
  filter(adj.P.Val < FDR_threshold)

deg_up <- deg_significant %>%
  filter(logFC > 0)

deg_down <- deg_significant %>%
  filter(logFC < 0)

# -----------------------------------------------------
# Save DEG subsets
# -----------------------------------------------------

write.csv(
  deg_significant,
  "results/DEG_significant_AD_vs_control.csv",
  row.names = FALSE
)

write.csv(
  deg_up,
  "results/DEG_upregulated_AD_vs_control.csv",
  row.names = FALSE
)

write.csv(
  deg_down,
  "results/DEG_downregulated_AD_vs_control.csv",
  row.names = FALSE
)

# -----------------------------------------------------
# Final integrity checks
# -----------------------------------------------------

cat("\nVerifying differential expression output integrity...\n")

if (!all(c("logFC", "adj.P.Val", "GeneSymbol") %in% colnames(results))) {
  stop("Expected columns missing in limma results.")
}

if (nrow(results) == 0) {
  stop("No genes detected in differential expression output.")
}

cat("Output structure verified.\n")

# -----------------------------------------------------
# Summary snapshot
# -----------------------------------------------------

cat("\n============================================\n")
cat("Differential Expression Analysis Complete\n")
cat("============================================\n")
cat("Total genes tested:", nrow(results), "\n")
cat("Significant genes (FDR < 0.05):", nrow(deg_significant), "\n")
cat("Upregulated:", nrow(deg_up), "\n")
cat("Downregulated:", nrow(deg_down), "\n")
cat("Results saved to: results/\n")
cat("============================================\n")

# -----------------------------------------------------
# End of Script
# -----------------------------------------------------

cat("\nScript 02 finished successfully.\n")

# =====================================================
# END OF SCRIPT 02
# Differential expression layer established
# Transcriptomic axis now defined for downstream analysis
# =====================================================


