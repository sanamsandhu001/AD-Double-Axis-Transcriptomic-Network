


# ============================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 00_setup_packages.R
# Author: Sanampreet Kaur
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 600)
set.seed(123)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/enrichment", showWarnings = FALSE)
dir.create("results/module_enrichment", showWarnings = FALSE)

cran_packages <- c(
  "dplyr", "ggplot2", "pheatmap",
  "igraph", "readr", "stringr"
)

bioc_packages <- c(
  "GEOquery", "limma",
  "clusterProfiler", "org.Hs.eg.db", "STRINGdb"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(igraph)
library(readr)
library(stringr)

cat("Environment ready.\n")

# ============================================================
# END OF SCRIPT 00
# Environment ready. Reproducible analysis ensured.
# ============================================================







