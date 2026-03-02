


# =====================================================
# Project: Double-Axis Transcriptomic Model of Alzheimer’s Disease
# Script: 08_publication_figures.R
# Author: Sanampreet Kaur
# Description: Refined publication-grade figures
# =====================================================

library(ggplot2)
library(dplyr)
library(igraph)

dir.create("figures", showWarnings = FALSE)

# =====================================================
# FIGURE 1 — Refined Volcano Plot
# =====================================================

limma_path <- "results/limma_full_results_AD_vs_control.csv"

if (!file.exists(limma_path)) {
  stop("Limma results not found. Run Script 02 first.")
}

limma_results <- read.csv(limma_path, stringsAsFactors = FALSE)

limma_results$Significance <- "NS"
limma_results$Significance[limma_results$adj.P.Val < 0.05 &
                             limma_results$logFC > 1] <- "Up"
limma_results$Significance[limma_results$adj.P.Val < 0.05 &
                             limma_results$logFC < -1] <- "Down"

volcano_plot <- ggplot(limma_results,
                       aes(x = logFC,
                           y = -log10(adj.P.Val))) +
  geom_point(data = subset(limma_results, Significance == "NS"),
             color = "grey85", alpha = 0.4, size = 1) +
  geom_point(data = subset(limma_results, Significance == "Up"),
             color = "#D55E00", alpha = 0.8, size = 1.2) +
  geom_point(data = subset(limma_results, Significance == "Down"),
             color = "#0072B2", alpha = 0.8, size = 1.2) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 15) +
  labs(title = "Transcriptomic Remodeling in AD Cortex",
       x = "Log2 Fold Change",
       y = "-log10 FDR")

ggsave("figures/Figure1_Volcano.pdf",
       volcano_plot,
       width = 6.5, height = 5.5)

# =====================================================
# FIGURE 2 — Refined Downregulated Network
# =====================================================

edges_path <- "data/PPI_Downregulated_interactions.csv"
modules_path <- "results/Louvain_Downregulated_Module_Assignments.csv"

if (!file.exists(edges_path) | !file.exists(modules_path)) {
  stop("Required PPI or Louvain outputs not found. Run Scripts 05 and 06 first.")
}

edges_down <- read.csv(edges_path, stringsAsFactors = FALSE)
modules_down <- read.csv(modules_path, stringsAsFactors = FALSE)

# Detect correct edge column names automatically
if (all(c("protein1_symbol", "protein2_symbol") %in% colnames(edges_down))) {
  colnames(edges_down)[colnames(edges_down) == "protein1_symbol"] <- "from"
  colnames(edges_down)[colnames(edges_down) == "protein2_symbol"] <- "to"
}

top_modules <- modules_down %>%
  count(Module, sort = TRUE) %>%
  slice(1:3) %>%
  pull(Module)

selected_nodes <- modules_down %>%
  filter(Module %in% top_modules)

edges_filtered <- edges_down %>%
  filter(from %in% selected_nodes$GeneSymbol &
           to %in% selected_nodes$GeneSymbol)

g_down <- graph_from_data_frame(edges_filtered, directed = FALSE)

layout_coords <- layout_with_fr(g_down)

module_vector <- selected_nodes$Module[
  match(V(g_down)$name, selected_nodes$GeneSymbol)
]

module_factor <- factor(module_vector)
module_colors <- c("#0072B2", "#D55E00", "#009E73")
color_map <- setNames(module_colors[1:length(levels(module_factor))],
                      levels(module_factor))

pdf("figures/Figure2_Downregulated_Network.pdf",
    width = 7, height = 6)

plot(g_down,
     layout = layout_coords,
     vertex.size = 6,
     vertex.label = NA,
     vertex.color = color_map[as.character(module_factor)],
     edge.color = adjustcolor("grey70", alpha.f = 0.4),
     main = "Downregulated Network — Top 3 Modules")

legend("topright",
       legend = paste("Module", levels(module_factor)),
       col = color_map,
       pch = 19,
       pt.cex = 1.5,
       bty = "n")

dev.off()

# =====================================================
# FIGURE 3 — Refined Modularity Comparison
# =====================================================

robust_down_path <- "results/Louvain_Downregulated_Robustness_Stats.csv"
robust_up_path   <- "results/Louvain_Upregulated_Robustness_Stats.csv"

if (!file.exists(robust_down_path) | !file.exists(robust_up_path)) {
  stop("Robustness stats not found. Run Script 06 first.")
}

robust_down <- read.csv(robust_down_path)
robust_up   <- read.csv(robust_up_path)

mod_df <- data.frame(
  Network = c("Downregulated", "Downregulated",
              "Upregulated", "Upregulated"),
  Type = c("Real", "Random", "Real", "Random"),
  Q = c(robust_down$Real_Modularity,
        robust_down$Random_Modularity,
        robust_up$Real_Modularity,
        robust_up$Random_Modularity)
)

mod_plot <- ggplot(mod_df,
                   aes(x = Network,
                       y = Q,
                       fill = Type)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.6),
           width = 0.6) +
  geom_text(aes(label = round(Q, 2)),
            position = position_dodge(width = 0.6),
            vjust = -0.4,
            size = 4) +
  scale_fill_manual(values = c("Real" = "#0072B2",
                               "Random" = "grey60")) +
  ylim(0, max(mod_df$Q) + 0.1) +
  theme_minimal(base_size = 15) +
  labs(title = "Network Modularity vs Random",
       y = "Modularity (Q)",
       x = "")

ggsave("figures/Figure3_Modularity_Comparison.pdf",
       mod_plot,
       width = 6.5, height = 5.5)

# =====================================================
# END OF SCRIPT 08
# Refined publication figures generated
# =====================================================