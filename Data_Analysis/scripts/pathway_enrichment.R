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


filtered_data_file <- "data/de_results_filtered.csv"

# Read in filtered DEGs
filtered <- read.csv(filtered_data_file, row.names = 1)
significant_genes <- rownames(filtered)

# GO Enrichment
go_results <- enrichGO(
  gene          = significant_genes, 
  OrgDb         = org.Hs.eg.db, 
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH", 
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
#write.csv(as.data.frame(go_results), "data/go_results.csv", row.names = FALSE)

# Visualize GO Enrichment
barplot(go_results, showCategory = 10, title = "Top 10 Enriched GO Terms (Biological Process)")

dotplot(go_results, showCategory = 10) +
  ggtitle("Dotplot of GO Enrichment (Biological Process)")

# Network plot for top 5 enriched terms
cnetplot(go_results, showCategory = 5) +
  ggtitle("GO Enrichment Network Plot")

# Combination plot
go_bar <- barplot(go_results, showCategory = 10, title = "Top 10 Enriched GO Terms (Biological Process)")
go_dot <- dotplot(go_results, showCategory = 10) + ggtitle("Dotplot of GO Enrichment (Biological Process)")
go_network <- cnetplot(go_results, showCategory = 5) + 
  ggtitle("GO Enrichment Network Plot") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- grid.arrange(
  go_bar, go_dot, 
  go_network, 
  ncol = 2, nrow = 2,
  heights = c(1, 2),
  layout_matrix = matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
)
print(combined_plot)

# KEGG Pathway Enrichment Analysis
# Convert SYMBOL to ENTREZID
gene_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG enrichment
kegg_results <- enrichKEGG(
  gene         = gene_ids$ENTREZID, 
  organism     = 'hsa',
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
#write.csv(as.data.frame(kegg_results), "data/kegg_results.csv", row.names = FALSE)

# Visualize KEGG Enrichment
barplot(kegg_results, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")

dotplot(kegg_results, showCategory = 10) +
  ggtitle("Dotplot of KEGG Enrichment")

# Reactome Pathway Analysis
reactome_results <- enrichPathway(
  gene = gene_ids$ENTREZID, 
  organism = "human", 
  pvalueCutoff = 0.05
)
#write.csv(as.data.frame(reactome_results), "data/reactome_results.csv", row.names = FALSE)

# Visualize Reactome Pathway Enrichment
barplot(reactome_results, showCategory = 10, title = "Top 10 Enriched Reactome Pathways")

dotplot(reactome_results, showCategory = 10) +
  ggtitle("Dotplot of Reactome Enrichment")
