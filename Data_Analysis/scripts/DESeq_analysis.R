#### DATA PREPARATION ####
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

#### Load all required libraries and ensure data ready for DESeq ####
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(apeglm)
library(BiocManager)
library(patchwork)
library(cowplot)

# Load count data
counts<-read.csv('data/GSE217504_host_counts_matrix.csv', header = T,row.names = 1)
colnames(counts)

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

#### PERFORMING DESEQ  ####
#create deseq object (this produces a warning about dropping factor levels - this refers to the hours no longer in use, we're only looking at hours 4, 12 and 48 as that is all we have mock data for)
dds<- DESeqDataSetFromMatrix(countData = count_data, colData = samples, design = ~Hours + Treatment)

# Set the reference for the Treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c("mock", "infected"))

# Filter the genes to only preserve those with 5 counts or more
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

# Extract the most differetially expresed genes due to the Treatment.
# select genes with a significant change in gene expression (adjusted p-value below 0.05)
# And log2fold change <-1 and >1
filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)

# Save the deseq reults
#write.csv(deseq_result,'data/de_results_all.csv')
#write.csv(filtered,'data/de_results_filtered.csv')

# Save the normalised counts
normalised_counts <- counts(dds,normalized=T)
#write.csv(normalised_counts,'data/normalised_counts.csv')

#### Repeat DESEQ for every timepoint
samples_4h <- subset(samples, Hours == "4")
count_data_4h <- count_data[, rownames(samples_4h)]
dds_4h <- DESeqDataSetFromMatrix(countData = count_data_4h, colData = samples_4h, design = ~ Treatment)
dds_4h$Treatment <- factor(dds_4h$Treatment, levels = c("mock", "infected"))
keep_4h <- rowSums(counts(dds_4h)) >= 5
dds_4h <- dds_4h[keep,]
dds_4h <- DESeq(dds_4h)
res_4h <- results(dds_4h)
res_4h <- as.data.frame(res_4h)
#write.csv(res_4h, 'data/de_results_4h.csv')
normalised_counts_4h <- counts(dds_4h,normalized = T)
#write.csv(normalised_counts_4h, 'data/normalised_counts_4h.csv')
filtered_4h <- res_4h %>% filter(res_4h$padj < 0.05)
filtered_4h <- filtered_4h %>% filter(abs(filtered_4h$log2FoldChange) > 1)
#write.csv(filtered_4h,'data/de_results_4h_filtered.csv')

samples_12h <- subset(samples, Hours == "12")
count_data_12h <- count_data[, rownames(samples_12h)]
dds_12h <- DESeqDataSetFromMatrix(countData = count_data_12h, colData = samples_12h, design = ~ Treatment)
dds_12h$Treatment <- factor(dds_12h$Treatment, levels = c("mock", "infected"))
keep_12h <- rowSums(counts(dds_12h)) >= 5
dds_12h <- dds_12h[keep,]
dds_12h <- DESeq(dds_12h)
res_12h <- results(dds_12h)
res_12h <- as.data.frame(res_12h)
#write.csv(res_12h, 'data/de_results_12h.csv')
normalised_counts_12h <- counts(dds_12h,normalized = T)
#write.csv(normalised_counts_12h, 'data/normalised_counts_12h.csv')
filtered_12h <- res_12h %>% filter(res_12h$padj < 0.05)
filtered_12h <- filtered_12h %>% filter(abs(filtered_12h$log2FoldChange) > 1)
#write.csv(filtered_12h,'data/de_results_12h_filtered.csv')

samples_48h <- subset(samples, Hours == "48")
count_data_48h <- count_data[, rownames(samples_48h)]
dds_48h <- DESeqDataSetFromMatrix(countData = count_data_48h, colData = samples_48h, design = ~ Treatment)
dds_48h$Treatment <- factor(dds_48h$Treatment, levels = c("mock", "infected"))
keep_48h <- rowSums(counts(dds_48h)) >= 5
dds_48h <- dds_48h[keep,]
dds_48h <- DESeq(dds_48h)
res_48h <- results(dds_48h)
res_48h <- as.data.frame(res_48h)
#write.csv(res_48h, 'data/de_results_48h.csv')
normalised_counts_48h <- counts(dds_48h,normalized = T)
#write.csv(normalised_counts_48h, 'data/normalised_counts_48h.csv')
filtered_48h <- res_48h %>% filter(res_48h$padj < 0.05)
filtered_48h <- filtered_48h %>% filter(abs(filtered_48h$log2FoldChange) > 1)
#write.csv(filtered_48h,'data/de_results_48h_filtered.csv')

#### EXPLORING & VISUALISING THE DATA ####

# Replace rownames as descriptive (I stopped using this as I liked the look of Hours/Treatment annotation for the heatmaps)
samples$SampleNumber <- ave(1:nrow(samples),samples$Treatment,samples$Hours,FUN = seq_along)
samples$SampleID <- paste0(samples$Treatment, "_",samples$Hours, "h_",samples$SampleNumber)
samples$SampleNumber <- NULL
print(samples)

# Dispersion plot
plotDispEsts(dds,main="Dispersion Estimates of Gene Expression")

# PCA plot
# variance stabalising transformation
vsd <- vst(dds,blind=F)

#use transformed values to generate a pca plot
pca_plot <- plotPCA(vsd, intgroup = c("Hours","Treatment"))
pca_plot + ggtitle("Principal Component Analysis (PCA) of Samples")

# Heatmap
#generate distance martrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <-as.matrix(sampleDists)

#set a colour scheme
colours <- colorRampPalette(rev(brewer.pal(9,"Greens")))(255)

annot_info <- as.data.frame(colData(dds)[,c('Treatment','Hours')])
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  color = colours,
  annotation_col = annot_info,
  main = "Heatmap of Sample-to-Sample Distances",
  #labels_col = samples$SampleID,
  #labels_row = samples$SampleID
)

# Heatmap of log transformed, using top 10 genes
top_hits <- deseq_result[order(deseq_result$padj),][1:10,]
top_hits <- row.names(top_hits)
top_hits

rld <- rlog(dds,blind=F)

label_colors <- ifelse(dds$Treatment == "mock", "blue", "red")

pheatmap(
  assay(rld)[top_hits,],
  annotation_col = annot_info,
  main = "Heatmap of Top 10 Most Expressed Genes",
  #labels_col = samples$SampleID
)

# Heatmap of Z scores. using top 10 genes.
cal_z_score <- function(x) {(x-mean(x))/sd(x)}

zscore_all <- t(apply(normalised_counts,1,cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset,
         annotation_col = annot_info,
         main = "Heatmap of Z-Scored Expression Levels for Top 10 Genes")

# Hierarchical Clustering Dendrogram
hc <- hclust(sampleDists, method = "ward.D2")
plot(hc, main = "Hierarchical Clustering Dendrogram of Samples", xlab = "", sub = "", cex = 0.8)
plot(hc, 
     labels = samples$SampleID,
     main = "Hierarchical Clustering Dendrogram of Samples",
     xlab = "Samples (Hours | Treatment)",
     sub = "",
     hang = -1)
rect.hclust(hc, k=3)

# MA Plot
plotMA(dds,ylim=c(-2,2))

#remove noise
resultsNames(dds)
resLFC <- lfcShrink(dds,coef="Treatment_infected_vs_mock", type="apeglm")

plotMA(resLFC,ylim=c(-2,2),main="MA Plot of Differential Gene Expression")

# Volcano Plots
resLFC <- as.data.frame(resLFC)

#label genes based on differential gene expression
resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05] <- "DOWN"
significant_genes <- resLFC[abs(resLFC$log2FoldChange) > 2 & resLFC$padj < 0.05, ]

# Label significant genes
resLFC$delabel <- NA
resLFC$delabel[rownames(resLFC) %in% rownames(significant_genes)] <- rownames(significant_genes)

# Volcano plot
ggplot(data=resLFC, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel)) +
  geom_point(aes(size = -log10(padj)), alpha=0.8) +
  theme_minimal() +
  geom_text_repel(data=subset(resLFC, abs(log2FoldChange) > 2 & padj < 0.05),max.overlaps = 10,color="darkgreen") +
  scale_color_manual(values=c('skyblue', 'grey80', 'salmon'), name="Expression Change") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change (Infected vs Mock)", y = "-log10(adjusted p-value)",
       caption = "Threshold: padj < 0.05, Log2 Fold Change > |1|") +
  theme(text = element_text(size = 16), legend.position = "bottom")

# Created function to do the above volcano plot on all 3 time points
LFC_4h <- as.data.frame(lfcShrink(dds_4h,coef="Treatment_infected_vs_mock", type="apeglm"))
LFC_12h <- as.data.frame(lfcShrink(dds_12h,coef="Treatment_infected_vs_mock", type="apeglm"))
LFC_48h <- as.data.frame(lfcShrink(dds_48h,coef="Treatment_infected_vs_mock", type="apeglm"))

generate_volcano <- function(data, title) {
  
  data <- data[!is.na(data$log2FoldChange) & !is.na(data$padj), ]
  
  data$diffexpressed <- "NO"
  data$diffexpressed[data$log2FoldChange > 1 & data$padj < 0.05] <- "UP"
  data$diffexpressed[data$log2FoldChange < -1 & data$padj < 0.05] <- "DOWN"
  
  significant_genes <- data[abs(data$log2FoldChange) > 2 & data$padj < 0.05, ]
  data$delabel <- NA
  data$delabel[rownames(data) %in% rownames(significant_genes)] <- rownames(significant_genes)
  
  p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
    geom_point(aes(size = -log10(padj)), alpha = 0.8) +
    geom_text_repel(data = subset(data, abs(log2FoldChange) > 2 & padj < 0.05), max.overlaps = 10, color = "darkgreen") +
    scale_color_manual(values = c('skyblue', 'grey80', 'salmon'), name = "Expression Change") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
    labs(title = title, x = "Log2 Fold Change (Infected vs Mock)", y = "-log10(adjusted p-value)") +
    theme_minimal() +
    theme(text = element_text(size = 16), 
          legend.position = "none",
          axis.text = element_text(size = 12))
  
  return(p)
}

# Generate volcano plots for 4h, 12h, and 48h
plot_4h <- generate_volcano(LFC_4h, "Volcano Plot for 4h")
plot_12h <- generate_volcano(LFC_12h, "Volcano Plot for 12h")
plot_48h <- generate_volcano(LFC_48h, "Volcano Plot for 48h")

# Combine plots side by side
combined_plot <- (plot_4h + plot_12h + plot_48h) +
  plot_layout(ncol = 3, guides = 'collect') +
  plot_annotation(
    title = "Differential Gene Expression Over Time of Infected cells compared to Mock",
    subtitle = "Threshold: padj < 0.05, Log2 Fold Change > |1|",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  ) + 
  theme(legend.position = "bottom",legend.box.spacing = unit(1, "lines")) 

# Display the combined plot
combined_plot