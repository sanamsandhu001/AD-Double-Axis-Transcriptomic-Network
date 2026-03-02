



# ============================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 01_data_preparation.R
# Author: Sanampreet Kaur
# Description: Download, clean, normalize and collapse probes
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)

dir.create("data", showWarnings = FALSE)

# ------------------------------------------------------------
# Download GEO dataset
# ------------------------------------------------------------

gset <- getGEO("GSE5281", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

saveRDS(gset, "data/gset_GSE5281.rds")

# ------------------------------------------------------------
# Extract expression matrix
# ------------------------------------------------------------

expr_matrix <- exprs(gset)

# ------------------------------------------------------------
# Log2 transform if needed
# ------------------------------------------------------------

qx <- quantile(expr_matrix, c(0, 0.25, 0.5, 0.75, 1))

if (qx[5] > 100) {
  expr_matrix <- log2(expr_matrix + 1)
}

expr_matrix <- normalizeBetweenArrays(expr_matrix)

# ------------------------------------------------------------
# Extract annotation
# ------------------------------------------------------------

annotation <- fData(gset)

symbol_column <- grep("symbol",
                      colnames(annotation),
                      ignore.case = TRUE,
                      value = TRUE)[1]

annotation$GeneSymbol <- annotation[[symbol_column]]

# Remove probes without gene symbol
valid_probes <- !is.na(annotation$GeneSymbol) &
  annotation$GeneSymbol != ""

expr_matrix <- expr_matrix[valid_probes, ]
annotation  <- annotation[valid_probes, ]

# ------------------------------------------------------------
# Collapse multiple probes per gene (FIXED SECTION)
# ------------------------------------------------------------

expr_df <- as.data.frame(expr_matrix)
expr_df$GeneSymbol <- annotation$GeneSymbol

expr_collapsed <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Convert to proper matrix and preserve rownames
expr_collapsed_mat <- as.matrix(expr_collapsed[, -1])
rownames(expr_collapsed_mat) <- expr_collapsed$GeneSymbol

cat("Unique genes after collapsing:", nrow(expr_collapsed_mat), "\n")

# ------------------------------------------------------------
# Extract phenotype data
# ------------------------------------------------------------

pdata <- pData(gset)

group_column <- "disease state:ch1"

pdata$Group <- ifelse(
  grepl("Alzheimer", pdata[[group_column]], ignore.case = TRUE),
  "AD",
  "Control"
)

pdata$Group <- factor(pdata$Group)

# Align samples
common_samples <- intersect(colnames(expr_collapsed_mat),
                            rownames(pdata))

expr_collapsed_mat <- expr_collapsed_mat[, common_samples]
pdata <- pdata[common_samples, ]

# ------------------------------------------------------------
# Save cleaned data
# ------------------------------------------------------------

saveRDS(expr_collapsed_mat, "data/expr_GSE5281_clean.rds")
saveRDS(pdata, "data/pdata_GSE5281_clean.rds")

cat("Cleaned gene-level expression saved.\n")

# ============================================================
# END OF SCRIPT 01
# Gene-level normalized dataset ready
# ============================================================

