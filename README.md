# AD-Double-Axis-Transcriptomic-Network

[![GitHub release](https://img.shields.io/github/v/release/sanamsandhu001/AD-Double-Axis-Transcriptomic-Network)](https://github.com/sanamsandhu001/AD-Double-Axis-Transcriptomic-Network/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Integrative Transcriptomic and Protein–Protein Interaction Network Analysis Reveals Divergent Mitochondrial–Synaptic Suppression and Immune Network Activation in Alzheimer's Disease**

This repository contains the R and Python analysis pipelines, network models, centrality data, and functional enrichment workflows supporting the systems-level dissection of Alzheimer's disease (AD) entorhinal cortex pathology.

---

## Overview

Bulk transcriptomic profiling of AD post-mortem tissue routinely detects thousands of differentially expressed genes, making it challenging to differentiate upstream drivers from reactive consequences. This pipeline integrates:

1. **Differential Expression Analysis:** Linear modeling via `limma` on laser-capture microdissected entorhinal cortex neurons (GEO: [GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281)).
2. **Protein–Protein Interaction (PPI) Networks:** Independent high-confidence interactomes (STRING v11.5, combined score ≥ 0.700) built for downregulated and upregulated axes.
3. **Unsupervised Community Detection:** Louvain modularity optimization benchmarked across repeated stochastic runs and matched Erdős–Rényi random graphs ($G(n, m)$).
4. **Resampling-Validated Hub Prioritization:** Dual-centrality candidate filtering (top decile of degree and normalized betweenness) tested against a 100-iteration, 10% edge-dropout bootstrap model.
5. **GWAS Integration:** Systematic intersection of network hubs against genome-wide significant ($P < 5 \times 10^{-8}$) loci from the NHGRI-EBI GWAS Catalog.

---

## Key Findings

* **Distinct Transcriptomic Scales:** 4,215 genes reached significance ($\text{adj. } P < 0.05$); downregulated genes (2,585 genes; 61.3%) significantly outnumbered upregulated genes (1,630 genes; 38.7%).
* **Topological Non-Randomness:** Louvain clustering identified 20 communities in the downregulated interactome ($Q = 0.751$ vs. $0.351$ random null) and 21 communities in the upregulated interactome ($Q = 0.709$ vs. $0.439$ random null).
* **Functional Modularity:**
  * **Downregulated:** Synaptic vesicle cycle & transmission (M6, M9), ubiquitin–proteasome catabolism (M2), DNA repair & transcription (M7), and GPI/glycan biosynthesis (M5).
  * **Upregulated:** Epigenetic regulation (M1), cytokine/JAK-STAT/PI3K-Akt signaling (M9), ECM/integrin signaling (M8), and innate immunity via TLR/NF-κB (M4).
* **Hub Prioritization & Genetic Triangulation:** Identified 55 robust downregulated hubs (e.g., `NDUFB7`, `ATP5A1`, `ATP5B`, `SNAP25`, `GRIN1`) and 23 robust upregulated hubs (e.g., `EP300`, `NFKB1`, `IL6`, `PIK3R1`).
* **GWAS Divergence:** Only a single robust hub gene—**`ATP5C1`** (downregulated mitochondrial ATP synthase subunit, bootstrap stability 0.98)—converged directly with AD GWAS-implicated risk loci, demonstrating that network hubs predominantly represent convergent downstream disease programs rather than primary germline risk drivers.

---

## Repository Structure

```text
├── data/
│   ├── raw/                           # Raw expression & metadata downloads
│   └── processed/                     # Quantile-normalized expression matrices
├── scripts/
│   ├── 01_data_preprocessing.R        # GEO download, probe summarization, normalization
│   ├── 02_differential_expression.R   # limma empirical Bayes linear modeling
│   ├── 03_ppi_network_construction.R  # STRINGdb mapping & igraph component extraction
│   ├── 04_louvain_modularity.R        # Community detection & Erdos-Renyi null comparisons
│   ├── 05_functional_enrichment.R     # clusterProfiler GO BP & KEGG over-representation
│   ├── 06_bootstrap_hub_validation.py # NetworkX bootstrap edge-dropout (100 iterations)
│   └── 07_gwas_intersection.R         # NHGRI-EBI catalog curation & hub matching
├── figures/                           # Generated publication figures (PDF / PNG)
│   ├── Figure1A_Volcano_Plot.png
│   ├── Figure1B_Top50_Heatmap.png
│   ├── Figure3_Modularity_Null_Plot.png
│   ├── Figure5_Downregulated_Interactome.png
│   ├── Figure6_Hub_Centrality_GWAS.png
│   └── Figure7_Module_Enrichment_Dotplot.png
├── results/                           # Output metric tables, gene lists, and statistics
│   ├── Downregulated_Hub_GWAS_Overlap.csv
│   ├── Upregulated_Hub_GWAS_Overlap.csv
│   └── PPI_All_Nodes_Centrality_GWAS.csv
├── GWAS_AD_associations.tsv           # Curated genome-wide significant AD associations
├── Transcriptomic Analysis of AD.Rproj# RStudio Project environment
├── .gitignore                         # Local caches & runtime artifacts
└── README.md
## Citation

If you use this codebase, network datasets, or prioritized hub targets in your research, please cite both the software archive and the accompanying study:

### Software Archive (Zenodo)
```bibtex
@software{kaur_2026_zenodo_23529818,
  author       = {Kaur, Sanampreet},
  title        = {AD-Double-Axis-Transcriptomic-Network: A double axis transcriptomic network model of Alzheimer's disease cortex},
  month        = oct,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v2.1},
  doi          = {10.5281/zenodo.23529818},
  url          = {[https://doi.org/10.5281/zenodo.23529818](https://doi.org/10.5281/zenodo.23529818)}
}
