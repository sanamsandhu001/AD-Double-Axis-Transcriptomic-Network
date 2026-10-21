# ==============================================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script:  08_Publication_Figures_Pipeline.R
# Author:  Sanampreet Kaur
# Description: End-to-end publication figure generation pipeline
#              - Figure 1: Volcano Plot (Transcriptomic Remodeling)
#              - Figure 2: Top 50 DEG Hierarchical Heatmap
#              - Figure 3: Downregulated PPI Network Core Modules
#              - Figure 4: Topological Modularity Validation (Real vs. Null)
#              - Figure 5: Functional Pathway Enrichment Dot Plot
#
# NOTE: all dplyr verbs below (filter, mutate, arrange, slice, slice_min,
# count, pull, recode, bind_rows, case_when) are explicitly namespaced as
# dplyr::verb(). Without this, slice() in particular resolves to the
# IRanges/S4Vectors S4 generic (defined for Rle objects) whenever those
# Bioconductor packages are attached in the same session -- e.g. via
# STRINGdb or clusterProfiler earlier in the pipeline -- producing:
#   "error in evaluating the argument 'x' in selecting a method for
#    function 'slice': Rle of type 'list' is not supported"
# This is the same collision class already fixed in scripts 03/05/07;
# it had not yet been applied here.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(igraph)
  library(readr)
  library(stringr)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
})

# Ensure output directory exists
dir.create("figures", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# SECTION 0: Upfront Input Validation
# ------------------------------------------------------------------------------
cat("\n=== [0/5] Checking Required Data Inputs ===\n")

required_files <- c(
  limma       = "results/limma_full_results_AD_vs_control.csv",
  expr        = "data/expr_GSE5281_clean.rds",
  pdata       = "data/pdata_GSE5281_clean.rds",
  edges_down  = "data/PPI_Downregulated_interactions.csv",
  modules_down= "results/Louvain_Downregulated_Module_Assignments.csv",
  robust_down = "results/Louvain_Downregulated_Robustness_Stats.csv",
  robust_up   = "results/Louvain_Upregulated_Robustness_Stats.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(paste0(
    "Missing required input files:\n",
    paste(paste0("  - ", names(missing_files), ": ", missing_files), collapse = "\n"),
    "\nPlease execute upstream pipeline scripts (01-07) before running Script 08."
  ))
}
cat("All core input dependencies verified.\n")


# ------------------------------------------------------------------------------
# SECTION 1: FIGURE 1 — Refined Volcano Plot
# ------------------------------------------------------------------------------
cat("\n=== [1/5] Generating Figure 1: Volcano Plot ===\n")

limma_results <- read.csv(required_files["limma"], stringsAsFactors = FALSE)

limma_results <- limma_results %>%
  dplyr::mutate(
    Significance = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE                          ~ "NS"
    )
  )

volcano_plot <- ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(data = dplyr::filter(limma_results, Significance == "NS"),
             color = "grey85", alpha = 0.4, size = 1.0) +
  geom_point(data = dplyr::filter(limma_results, Significance == "Down"),
             color = "#0072B2", alpha = 0.8, size = 1.3) +
  geom_point(data = dplyr::filter(limma_results, Significance == "Up"),
             color = "#D55E00", alpha = 0.8, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Transcriptomic Remodeling in AD Cortex",
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"FDR")
  )

ggsave("figures/Figure1_Volcano.pdf", volcano_plot, width = 6.5, height = 5.5)
ggsave("figures/Figure1_Volcano.png", volcano_plot, width = 6.5, height = 5.5, dpi = 300)
cat(" -> Saved: Figure1_Volcano (PDF/PNG)\n")


# ------------------------------------------------------------------------------
# SECTION 2: FIGURE 2 — Enhanced Top 50 DEG Heatmap
# ------------------------------------------------------------------------------
cat("\n=== [2/5] Generating Figure 2: Top 50 DEG Heatmap ===\n")

expr_matrix <- readRDS(required_files["expr"])
pheno_data  <- readRDS(required_files["pdata"])

top_50_df <- limma_results %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice(1:50) %>%
  dplyr::mutate(
    CleanSymbol = ifelse(grepl("///", GeneSymbol), 
                         str_extract(GeneSymbol, "[^/]+$"), 
                         GeneSymbol)
  )

valid_genes <- intersect(top_50_df$GeneSymbol, rownames(expr_matrix))
heatmap_mat <- expr_matrix[valid_genes, ]
rownames(heatmap_mat) <- top_50_df$CleanSymbol[match(rownames(heatmap_mat), top_50_df$GeneSymbol)]

annotation_col <- data.frame(
  Diagnosis = factor(pheno_data$Group, levels = c("Control", "AD")),
  row.names = rownames(pheno_data)
)

annotation_row <- data.frame(
  Regulation = factor(
    ifelse(top_50_df$logFC[match(rownames(heatmap_mat), top_50_df$CleanSymbol)] > 0, 
           "Upregulated", "Downregulated"),
    levels = c("Upregulated", "Downregulated")
  ),
  row.names = rownames(heatmap_mat)
)

ann_colors <- list(
  Diagnosis  = c("Control" = "#2B83BA", "AD" = "#D7191C"),
  Regulation = c("Upregulated" = "#E66101", "Downregulated" = "#5E3C99")
)

palette_expr <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)

render_heatmap <- function(file_path) {
  pheatmap(
    heatmap_mat,
    scale = "row",
    color = palette_expr,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 7.5,
    fontsize_col = 8,
    treeheight_row = 28,
    treeheight_col = 22,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    cutree_cols = 2,
    cutree_rows = 2,
    border_color = "grey92",
    main = "Top 50 Differentially Expressed Genes in AD Entorhinal Cortex",
    filename = file_path,
    width = 8.0,
    height = 10.5
  )
}

render_heatmap("figures/Figure2_Heatmap_Top50.pdf")
render_heatmap("figures/Figure2_Heatmap_Top50.png")
cat(" -> Saved: Figure2_Heatmap_Top50 (PDF/PNG)\n")


# ------------------------------------------------------------------------------
# SECTION 3: FIGURE 3 — Downregulated PPI Sub-Network Architecture
# ------------------------------------------------------------------------------
cat("\n=== [3/5] Generating Figure 3: Downregulated PPI Network ===\n")

edges_down   <- read.csv(required_files["edges_down"], stringsAsFactors = FALSE)
modules_down <- read.csv(required_files["modules_down"], stringsAsFactors = FALSE)
hubs_path    <- "results/Downregulated_Robust_Hubs.csv"

# Normalize edge column names
if (all(c("protein1_symbol", "protein2_symbol") %in% colnames(edges_down))) {
  colnames(edges_down)[colnames(edges_down) == "protein1_symbol"] <- "from"
  colnames(edges_down)[colnames(edges_down) == "protein2_symbol"] <- "to"
}

top_3_modules <- modules_down %>%
  dplyr::count(Module, sort = TRUE) %>%
  dplyr::slice(1:3) %>%
  dplyr::pull(Module)

module_labels <- c(
  "2" = "Module 2: Ubiquitin-Proteasome Catabolism",
  "6" = "Module 6: mTOR & Synaptic Vesicle Transport",
  "9" = "Module 9: Presynaptic Vesicle Cycle"
)

selected_nodes <- modules_down %>%
  dplyr::filter(Module %in% top_3_modules) %>%
  dplyr::mutate(Module_Name = dplyr::recode(as.character(Module), !!!module_labels))

edges_filtered <- edges_down %>%
  dplyr::filter(from %in% selected_nodes$GeneSymbol & to %in% selected_nodes$GeneSymbol)

g_sub <- graph_from_data_frame(edges_filtered, directed = FALSE, vertices = selected_nodes)
g_sub <- igraph::simplify(g_sub, remove.multiple = TRUE, remove.loops = TRUE)
V(g_sub)$degree <- degree(g_sub)

top_hubs <- if (file.exists(hubs_path)) {
  read.csv(hubs_path, stringsAsFactors = FALSE)$Gene[1:12]
} else {
  names(sort(degree(g_sub), decreasing = TRUE))[1:10]
}

palette_modules <- c(
  "Module 2: Ubiquitin-Proteasome Catabolism"   = "#3182BD",
  "Module 6: mTOR & Synaptic Vesicle Transport" = "#E6550D",
  "Module 9: Presynaptic Vesicle Cycle"         = "#31A354"
)

node_colors <- palette_modules[selected_nodes$Module_Name]
node_sizes  <- rescale(V(g_sub)$degree, to = c(4, 9.5))
node_labels <- ifelse(V(g_sub)$name %in% top_hubs, V(g_sub)$name, NA)

set.seed(123)
layout_fr <- layout_with_fr(g_sub, niter = 1200)

plot_network_device <- function() {
  par(mar = c(3, 2, 3, 2), bg = "white")
  plot(
    g_sub,
    layout = layout_fr,
    vertex.size = node_sizes,
    vertex.color = node_colors,
    vertex.frame.color = "white",
    vertex.frame.width = 0.5,
    vertex.label = node_labels,
    vertex.label.cex = 0.85,
    vertex.label.font = 2,
    vertex.label.color = "#111111",
    vertex.label.dist = 0.9,
    edge.color = adjustcolor("grey70", alpha.f = 0.35),
    edge.width = 0.7,
    main = "Downregulated Interactome: Core Sub-Modules"
  )
  legend(
    "bottom",
    legend = names(palette_modules),
    col = palette_modules,
    pch = 19,
    pt.cex = 1.6,
    bty = "n",
    cex = 0.9,
    inset = c(0, -0.04),
    xpd = TRUE
  )
}

pdf("figures/Figure3_Downregulated_Network.pdf", width = 8.5, height = 8.0)
plot_network_device()
dev.off()

png("figures/Figure3_Downregulated_Network.png", width = 8.5, height = 8.0, units = "in", res = 300)
plot_network_device()
dev.off()

cat(" -> Saved: Figure3_Downregulated_Network (PDF/PNG)\n")


# ------------------------------------------------------------------------------
# SECTION 4: FIGURE 4 — Modularity Validation vs. Null Models
# ------------------------------------------------------------------------------
cat("\n=== [4/5] Generating Figure 4: Modularity Comparison ===\n")

robust_down <- read.csv(required_files["robust_down"])
robust_up   <- read.csv(required_files["robust_up"])

mod_df <- data.frame(
  Network = factor(c("Downregulated", "Downregulated", "Upregulated", "Upregulated"),
                   levels = c("Downregulated", "Upregulated")),
  Topology = factor(c("Empirical (Real)", "Null Model (Random)", 
                      "Empirical (Real)", "Null Model (Random)"),
                    levels = c("Empirical (Real)", "Null Model (Random)")),
  Q = c(robust_down$Real_Modularity, robust_down$Random_Modularity,
        robust_up$Real_Modularity, robust_up$Random_Modularity),
  SD = c(robust_down$SD_Modularity_10_Runs, 0,
         robust_up$SD_Modularity_10_Runs, 0)
)

fig4_modularity <- ggplot(mod_df, aes(x = Network, y = Q, fill = Topology)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55, color = "#2B2B2B", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, Q - SD), ymax = Q + SD),
                position = position_dodge(width = 0.65), width = 0.18, linewidth = 0.6, color = "#2B2B2B") +
  geom_text(aes(label = sprintf("%.2f", Q)),
            position = position_dodge(width = 0.65), vjust = -1.2, size = 4.2, fontface = "bold", color = "#1F1F1F") +
  scale_fill_manual(values = c("Empirical (Real)" = "#1F4E79", "Null Model (Random)" = "#B0BEC5")) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.2), expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10.5),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 11),
    axis.line = element_line(color = "grey30", linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey35", margin = margin(b = 10))
  ) +
  labs(
    title = "Topological Community Modularity (Q) vs. Null Models",
    subtitle = "Mean ± SD across 10 Louvain runs vs. Erdos-Renyi random networks (matched n, m)",
    x = "Transcriptomic Axis",
    y = "Louvain Modularity Index (Q)"
  )

ggsave("figures/Figure4_Modularity_Comparison.pdf", fig4_modularity, width = 6.2, height = 5.2)
ggsave("figures/Figure4_Modularity_Comparison.png", fig4_modularity, width = 6.2, height = 5.2, dpi = 300)
cat(" -> Saved: Figure4_Modularity_Comparison (PDF/PNG)\n")


# ------------------------------------------------------------------------------
# SECTION 5: FIGURE 5 — Module Functional Enrichment Dot Plot
# ------------------------------------------------------------------------------
cat("\n=== [5/5] Generating Figure 5: Functional Enrichment Dot Plot ===\n")

enrich_dir <- "results/module_enrichment"
enrich_files <- list.files(enrich_dir, pattern = "_(GO_BP|KEGG)\\.csv$", full.names = TRUE)

if (length(enrich_files) > 0) {
  enrich_list <- list()
  for (f in enrich_files) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      fname    <- basename(f)
      net_type <- ifelse(grepl("Downregulated", fname), "Downregulated", "Upregulated")
      ont_type <- ifelse(grepl("GO_BP", fname), "GO BP", "KEGG")
      mod_num  <- str_extract(fname, "Module_([0-9]+)")
      
      df <- df %>%
        dplyr::slice_min(order_by = p.adjust, n = 3, with_ties = FALSE) %>%
        dplyr::mutate(
          Network = net_type,
          Ontology = ont_type,
          Module = mod_num,
          Module_Label = paste0(net_type, " ", mod_num),
          GeneRatio_Num = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
        )
      enrich_list[[fname]] <- df
    }
  }
  
  all_enrich_df <- dplyr::bind_rows(enrich_list)
  
  dotplot_fig5 <- ggplot(all_enrich_df, aes(x = GeneRatio_Num, y = reorder(Description, GeneRatio_Num))) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradient(low = "#4575b4", high = "#d73027", name = expression(-Log[10]~"(adj.P)")) +
    facet_grid(Module_Label ~ Ontology, scales = "free_y", space = "free_y") +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 8.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Functional Pathway Enrichment Across Network Modules",
      x = "Gene Ratio",
      y = "Pathway / Biological Process",
      size = "Gene Count"
    )
  
  ggsave("figures/Figure5_Module_Enrichment_Dotplot.pdf", dotplot_fig5, width = 11, height = 9)
  ggsave("figures/Figure5_Module_Enrichment_Dotplot.png", dotplot_fig5, width = 11, height = 9, dpi = 300)
  cat(" -> Saved: Figure5_Module_Enrichment_Dotplot (PDF/PNG)\n")
} else {
  cat(" [!] Warning: No module enrichment files found in results/module_enrichment/. Skipping Figure 5.\n")
}

cat("\n=======================================================\n")
cat(" Script 08 Pipeline Completed Successfully.\n")
cat(" All publication figures are available in ./figures/\n")
cat("=======================================================\n")