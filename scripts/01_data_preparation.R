# ============================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer's Disease
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

options(download.file.method = "wininet")  # Windows-only fix for schannel TLS errors

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
# Collapse multiple probes per gene
# ------------------------------------------------------------

expr_df <- as.data.frame(expr_matrix)
expr_df$GeneSymbol <- annotation$GeneSymbol

expr_collapsed <- expr_df %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), mean), .groups = "drop")

# Convert to proper matrix and preserve rownames
expr_collapsed_mat <- as.matrix(expr_collapsed[, -1])
rownames(expr_collapsed_mat) <- expr_collapsed$GeneSymbol

cat("Unique genes after collapsing:", nrow(expr_collapsed_mat), "\n")

# ------------------------------------------------------------
# Extract phenotype data
# ------------------------------------------------------------

pdata <- pData(gset)

# ------------------------------------------------------------
# Restrict to entorhinal cortex only, and define Group
# ------------------------------------------------------------
# NOTE: pData()'s "region" AND "disease state:ch1" characteristics
# columns are both unreliable for this dataset, for the same reason
# — GSE5281 spans multiple regions/submission batches, and
# GEOquery's characteristics parsing fragments labels across several
# columns depending on field order, leaving most samples as NA (or
# mismatched) in any single column. Both the EC restriction and the
# AD/Control classification are therefore done directly from the
# confirmed GSM accession lists (cross-verified against the GEO
# sample characteristics sheet), rather than by parsing pData.

ec_controls <- c(
  "GSM119615", "GSM119616", "GSM119617", "GSM119618", "GSM119619",
  "GSM119620", "GSM119621", "GSM119622", "GSM119623", "GSM119624",
  "GSM119625", "GSM119626", "GSM119627"
)

ec_ad <- c(
  "GSM238763", "GSM238790", "GSM238791", "GSM238792", "GSM238793",
  "GSM238794", "GSM238795", "GSM238796", "GSM238797", "GSM238798"
)

ec_samples <- c(ec_controls, ec_ad)

cat("Entorhinal cortex samples found:", sum(rownames(pdata) %in% ec_samples), "(expected: 23)\n")
stopifnot(sum(rownames(pdata) %in% ec_samples) == 23)

expr_collapsed_mat <- expr_collapsed_mat[, colnames(expr_collapsed_mat) %in% ec_samples]
pdata <- pdata[rownames(pdata) %in% ec_samples, ]

pdata$Group <- ifelse(
  rownames(pdata) %in% ec_ad, "AD",
  ifelse(rownames(pdata) %in% ec_controls, "Control", NA)
)
pdata$Group <- factor(pdata$Group, levels = c("Control", "AD"))

cat("Group counts:\n")
print(table(pdata$Group))
stopifnot(sum(pdata$Group == "Control") == 13, sum(pdata$Group == "AD") == 10)

# ------------------------------------------------------------
# Align samples
# ------------------------------------------------------------

common_samples <- intersect(colnames(expr_collapsed_mat),
                            rownames(pdata))

expr_collapsed_mat <- expr_collapsed_mat[, common_samples]
pdata <- pdata[common_samples, ]

stopifnot(all(colnames(expr_collapsed_mat) == rownames(pdata)))
cat("Final aligned sample count:", ncol(expr_collapsed_mat), "\n")

# ------------------------------------------------------------
# Save cleaned data
# ------------------------------------------------------------

saveRDS(expr_collapsed_mat, "data/expr_GSE5281_clean.rds")
saveRDS(pdata, "data/pdata_GSE5281_clean.rds")

cat("Cleaned gene-level expression saved.\n")

# ============================================================
# END OF SCRIPT 01
# Gene-level normalized dataset ready (entorhinal cortex only)
# ============================================================