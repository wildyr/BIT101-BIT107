# For data management
install.packages('tidyverse')
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
# For visualisation
install.packages('pheatmap')
install.packages("DOSE")
install.packages("enrichplot")
install.packages("ggupset")

library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(DOSE)
library(enrichplot)
library(ggupset)
library(dplyr)

df_4h <- read.csv("data/de_results_4h_filtered.csv")
df_12h <- read.csv("data/de_results_12h_filtered.csv")
df_48h <- read.csv("data/de_results_48h_filtered.csv")

colnames(df_4h)[1] <- "gene_symbol"
genes_in_4h <- df_4h$gene_symbol

pwl2 <- read.gmt("data/gmt/c5.go.bp.v2024.1.Hs.symbols.gmt")
pwl2 <- pwl2[pwl2$gene %in% genes_in_4h,]
saveRDS(pwl2, "data/gmt/go.bp.rds")

#gmt_files <- list.files(path = "data/gmt" , pattern = '.gmt', full.names = TRUE)
#for (file in gmt_files){
#  #file <- gmt_files[1]
#  pwl2 <- read.gmt(file) 
#  pwl2 <- pwl2[pwl2$gene %in% genes_in_4h,]
#  filename <- paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
#  saveRDS(pwl2, filename)
#}
#rm(file,filename,gmt_files)

df_4h <- df_4h %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0 & padj < 0.05 ~ 'UP',
  log2FoldChange < 0 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO'
))

df_4h <- df_4h[df_4h$diffexpressed != 'NO',]

deg_results_4h <- split(df_4h, df_4h$diffexpressed)


#### Running ClusterProfiler
#Settings
## Run ClusterProfiler -----------------------------------------------

# Settings
name_of_comparison <- 'infectedvsmock' # for our filename
background_genes <- 'go.bp' # for our filename
bg_genes <- readRDS(paste0('data/gmt/go.bp.RDS')) # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways
filename <- paste0("data/", 'clusterProfiler/', name_of_comparison, '_', background_genes) # filename of our PEA results
str(df_4h)
res_4h <- lapply(names(deg_results_4h),
                 function(x) enricher(gene = deg_results_4h[[x]]$gene_symbol,
                                      TERM2GENE = bg_genes))
names(res_4h) <- names(deg_results_4h)

res_4h_df <- lapply(names(res_4h), function(x) rbind(res_4h[[x]]@result))
names(res_4h_df) <- names(res_4h)
res_4h_df <- do.call(rbind, res_4h_df)
head(res_4h_df)

res_4h_df <- res_4h_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexpressed = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_4h_df)))

res_4h_df_up <- res_4h_df %>% filter(diffexpressed == 'UP') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexpressed'))
rownames(res_4h_df_up) <- res_4h_df_up$ID
colnames(res_4h_df_up)
# Assuming your data frame is called 'df'
res_4h_df_up <- res_4h_df_up[, !(colnames(res_4h_df_up) %in% c("RichFactor", "FoldEnrichment", "zScore"))]


enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_4h_df_up,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = df_4h$gene_symbol,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)
p1<- barplot(enrichres, showCategory = 20) 

p1
#target_pws <- unique(res_4h_df)

