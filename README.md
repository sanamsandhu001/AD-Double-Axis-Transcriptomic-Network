[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18833207.svg)](https://doi.org/10.5281/zenodo.18833207)

# AD-Double-Axis-Transcriptomic-Network

## Overview

This repository contains a fully reproducible transcriptomic network analysis pipeline for Alzheimer’s disease cortex.

The study proposes a double-axis systems model in which:

- Downregulated genes form structured synaptic modules exhibiting modular collapse.
- Upregulated genes reflect regulatory and architectural expansion.
- Network modularity significantly exceeds random expectation.

The pipeline integrates:

- Differential expression (limma)
- STRING protein–protein interaction networks
- Louvain community detection
- Random network comparison
- GO and KEGG enrichment
- Robust hub identification
- GWAS overlap analysis

---

## Dataset

- GEO accession: GSE5281  
- Tissue: Human cortex  
- Groups: Alzheimer’s disease vs Control  
- Species: Homo sapiens  

---

## Pipeline Structure

00_setup_packages.R  
01_data_preparation.R  
02_differential_expression_limma.R  
03_visualisation_plots.R  
04_functional_enrichment_GO_KEGG.R  
05_PPI_network_and_hub_genes.R  
06_Louvain_module_detection.R  
07_Module_Functional_Enrichment.R  
08_publication_figures.R  
09_GWAS_AD_overlap.R  

---

## Reproducibility

To reproduce the full analysis, run scripts sequentially from 00 to 09 in a fresh R session.

Session information is provided in `session_info.txt`.

---

## Conceptual Contribution

This project moves beyond gene-level differential expression toward network-structured interpretation of Alzheimer’s disease transcriptomic remodeling.

The results support a model of coordinated synaptic modular suppression embedded within reorganized interaction architecture.

---

## Author

Sanampreet Kaur
