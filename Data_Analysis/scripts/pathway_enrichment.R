#### Pathway Enrichment Analysis ####
# Install and load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

necessary_libraries <- c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot", "ggplot2", "dplyr", "ggpubr", "gridExtra")
for (lib in necessary_libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    BiocManager::install(lib)
  }
}
rm(necessary_libraries)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)

# Read in filtered DEGs
filtered <- read.csv("data/de_results_filtered.csv", row.names = 1)
upregulated <- filtered %>% filter(padj < 0.05 & log2FoldChange > 1)
downregulated <- filtered %>% filter(padj < 0.05 & log2FoldChange < -1)
up_significant_genes <- rownames(upregulated)
down_significant_genes <- rownames(downregulated)

# GO Enrichment
go_results_down <- enrichGO(
  gene          = down_significant_genes, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
go_results_up <- enrichGO(
  gene          = up_significant_genes, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.2
)

#write.csv(as.data.frame(go_results), "data/go_results.csv", row.names = FALSE)

# Visualize GO Enrichment
barplot(go_results_up,showCategory=10,title="Top 10 Up-regulated GO Terms")
barplot(go_results_down,showCategory=10,title="Top 10 Down-regulated GO Terms")

dotplot(go_results_up, showCategory = 10)+ggtitle("Dotplot of GO Enrichment")
dotplot(go_results_down, showCategory = 10)+ggtitle("Dotplot of GO Enrichment")

# Network plot for top 5 enriched terms
cnetplot(go_results_up, showCategory = 5)+ggtitle("Up-regulated GO Enrichment Network Plot")+theme(plot.title = element_text(hjust = 0.5))

cnetplot(go_results_down, showCategory = 5)+ggtitle("Down-regulated GO Enrichment Network Plot")+theme(plot.title = element_text(hjust = 0.5))

# Combination plot
go_bar <- barplot(go_results_down, showCategory = 10, title = "Top 10 Down-regulated GO Terms ")
go_dot <- dotplot(go_results_down, showCategory = 10) + ggtitle("Dotplot of Down-regulated GO Enrichment")
go_network <- cnetplot(go_results_down, showCategory = 5) + 
  ggtitle("Down-regulated GO Enrichment Network Plot") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot_go <- grid.arrange(
  go_bar, go_dot, 
  go_network, 
  ncol = 2, nrow = 2,
  heights = c(1, 2),
  layout_matrix = matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
)