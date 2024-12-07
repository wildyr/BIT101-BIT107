metadata <- read.csv("data/final_metadata.csv", header = FALSE, stringsAsFactors = FALSE)

metadata <- read.csv("data/prepared_metadata2.csv", header = FALSE, stringsAsFactors = FALSE)

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

# Filter the infected data to keep only rows where Hours match the mock Hours (4, 12, or 48)
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

# Load count data
counts<-read.csv('data/GSE217504_host_counts_matrix.csv', header = T,row.names = 1)
colnames(counts)
head(counts)

colnames(samples)
samples

# Ensure that the samples dataframe is ordered by the row names (sample IDs)
samples <- samples[order(rownames(samples)), ]

# Subset counts to keep only the columns matching the sample IDs in samples
counts_filtered <- counts[, colnames(counts) %in% rownames(samples)]

# Reorder the columns in counts_filtered to match the order in samples
counts_filtered <- counts_filtered[, match(rownames(samples), colnames(counts_filtered))]

# Check if the columns in counts_filtered match the order in samples
all(colnames(counts_filtered) == rownames(samples))  # Should return TRUE

count_data<-counts_filtered
colnames(counts)
rownames(samples)

levels(samples$Hours)

rm(counts_filtered)

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
