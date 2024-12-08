#### Data Preparation ####
#Please setwd to Project directory
#setwd("~/R/git/assignment2/Data_Analysis")
metadata <- read.csv("data/final_metadata.csv", header = FALSE, stringsAsFactors = FALSE)

# Transpose the data
metadata_t <- as.data.frame(t(metadata))

# Assign the first row as column names
colnames(metadata_t) <- metadata_t[1, ]

# Remove the first row from the data
metadata_t <- metadata_t[-1, ]

# Move Sample description to position 1
metadata_t <- metadata_t[, c("Description", setdiff(colnames(metadata_t), "Description"))]

# Set row names to the Sample description column
rownames(metadata_t) <- metadata_t$`Description`

# Remove the Sample description column
metadata_t <- metadata_t[, -1]

# Set as factors
metadata_t$Treatment <- as.factor(metadata_t$Treatment)
metadata_t$Hours <- as.factor(metadata_t$Hours)

summary(metadata_t)
metadata<-metadata_t

# Separate the infected and mock rows
infected <- subset(metadata, Treatment == "infected")
mock <- subset(metadata, Treatment == "mock")

# Get the unique hours from the mock data (we know these are 4, 12, and 48 in your case)
mock_hours <- unique(mock$Hours)
mock_hours

# Filter the infected data to keep only rows where Hours match the mock Hours
filtered_infected <- infected[infected$Hours %in% mock_hours, ]

# Combine the filtered infected rows with all the mock rows
samples <- rbind(filtered_infected, mock)
rm(mock, infected, filtered_infected, mock_hours, metadata, metadata_t)

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
library(apeglm)


# Load count data
counts<-read.csv('data/GSE217504_host_counts_matrix.csv', header = T,row.names = 1)
colnames(counts)
head(counts)

samples
str(samples)

# Ensure that the samples dataframe is ordered by the row names (sample IDs)
samples <- samples[order(rownames(samples)), ]

# Subset counts to keep only the columns matching the sample IDs in samples
counts_filtered <- counts[, colnames(counts) %in% rownames(samples)]

# Reorder the columns in counts_filtered to match the order in samples
counts_filtered <- counts_filtered[, match(rownames(samples), colnames(counts_filtered))]

# Check if the columns in counts_filtered match the order in samples
all(colnames(counts_filtered) == rownames(samples))  # Should return TRUE

count_data<-counts_filtered
colnames(count_data)
rownames(samples)

rm(counts_filtered, counts)

#### DESEQ ####
#create deseq object (this produces a warning about dropping factor levels - this refers to the hours no longer in use, we're only looking at hours 4, 12 and 48 as that is all we have mock data for)
dds<- DESeqDataSetFromMatrix(countData = count_data, colData = samples, design = ~Hours + Treatment)

# Set the reference for the Treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c("mock", "infected"))

# Filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# Stats test to identify differentially expressed genes
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result

# Form results into dataframe
deseq_result <- as.data.frame(deseq_result)
head(deseq_result)

# Order result table by increasing p value
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
head(deseq_result_ordered)

# Some queries
# Is ZC3H12A gene differentially expressed?
deseq_result["ZC3H12A",]

# Extract the most differetially expresed genes due to the Treatment.
# select genes with a significant change in gene expression (adjusted p-value below 0.05)
# And log2fold change <1 and >1
filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

dim(deseq_result)
dim(filtered)

# Save the deseq reults. We will save both the original and the filtered one
write.csv(deseq_result,'data/de_results_all.csv')
write.csv(filtered,'data/de_results_filtered.csv')

# Save the normalised counts
normalised_counts <- counts(dds,normalized=T)
head(normalised_counts)
write.csv(normalised_counts,'data/normalised_counts.csv')

#### EXPLORING THE DATA ####
# Dispersion plot
plotDispEsts(dds)

# PCA plot
# variance stabalising transformation
vsd <- vst(dds,blind=F)

#use transformed values to generate a pca plot
plotPCA(vsd,intgroup=c("Hours", "Treatment"))

# Heatmap
#generate distance martrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <-as.matrix(sampleDists)
colnames(sampleDistMatrix)

#set a colour scheme
colours <- colorRampPalette(rev(brewer.pal(9,"Greens")))(255)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = colours,
  annotation_col = samples,
  main = "Distance Heatmap"
)

## Clearly highest similarity among the 48 hour bucket, regardless of treatment. Also between the mock treatment data, regardless of hours.

# Heatmap of log transformed, using top 10 genes
top_hits <- deseq_result[order(deseq_result$padj),][1:10,]
top_hits <- row.names(top_hits)
top_hits

rld <- rlog(dds,blind=F)

pheatmap(assay(rld)[top_hits,], cluster_rows=F,show_rownames=T,cluster_cols=F)
pheatmap(assay(rld)[top_hits,],)

annot_info <- as.data.frame(colData(dds)[,c('Hours','Treatment')])
pheatmap(assay(rld)[top_hits,],annotation_col = annot_info)


# Heatmap of Z scores. using top 10 genes.
cal_z_score <- function(x) {(x-mean(x))/sd(x)}

zscore_all <- t(apply(normalised_counts,1,cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset, annotation_col = annot_info)


# MA Plot
plotMA(dds,ylim=c(-2,2))

#remove noise
resultsNames(dds)
resLFC <- lfcShrink(dds,coef="Treatment_infected_vs_mock", type="apeglm")

plotMA(resLFC,ylim=c(-2,2))

# Volcano Plot
resLFC <- as.data.frame(resLFC)

#label genes
# Update thresholds for differential expression
resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05] <- "DOWN"

# Label significant genes
resLFC$delabel <- NA
resLFC$delabel[abs(resLFC$log2FoldChange) > 2 & resLFC$padj < 0.05] <- rownames(resLFC)

# Volcano plot
ggplot(data=resLFC, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point(aes(size = -log10(padj)), alpha=0.8) +
  theme_minimal() +
  geom_text_repel(data=subset(resLFC, abs(log2FoldChange) > 2 & padj < 0.05),max.overlaps = 10,color="darkgreen") +
  scale_color_manual(values=c('blue', 'grey80', 'red'), name="Expression Change") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change (Treated vs Mock)", y = "-log10(adjusted p-value)",
       caption = "Threshold: padj < 0.05, Log2 Fold Change > |1|") +
  theme(text = element_text(size = 16), legend.position = "bottom")



