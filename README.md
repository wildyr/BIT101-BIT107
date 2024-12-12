# assignment2
Assignment for BIT101 and BIT107 - Short communication style mini-project report

Data Retrieval:-

A host of scripts to filter the metadata file GSE217504_series_metadata_matrix.txt to the relevant data for analysis.

getVariables.sh - Retrieves all rownames and prints them to variable_list.txt
filterScript.sh - Looks in variable_list_edit.txt and filters our data to those specified - filtered_metadata.txt
NOTE: variable_list_edit.txt has been manually edited for this dataset. If running all scripts consecutively this file and the original metadata file are required.
convert.CSV.sh - Cleans row names up for better readability and converts the file to a csv - converted_metadata.csv
prepareForR.sh - Rewrites duplicate row names and removes superflous information from our data points - prepared_metadata.csv
cutDown.sh - Further filters applied to remove unecessary fields before importing to R for analysis - final_metadata.csv

For ease of use I've also added a pipline.sh to run this scripts consequtively. Please note that both the original metadata txt and the variable_list_edit.txt are required.

These scripts were written for the purpose of analysing this specific dataset, and would require significant changes for different datasets.

Data Analysis:-

DESeq_analysis.R
Data Preparation - Rearranges the data from the final_metadata.csv file we generated in Data Retrieval to be ready for DESeq. Also loads necessary libraries and reads in count data.
Performing DESeq - Performs DESeq analysis on the whole dataset and then separately on each different time point.
Exploring & Visualising the Data - Creates numerous plots on the DESeq of the whole dataset, followed by a volcano plot for each time point.

pathway_enrichment.R
Data from DESeq_analysis.R is used.

If running these scripts locally, please construct a file path to mirror to what you find in this git repository, or update all the file paths in the script.

Version Info:
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3               ggpubr_0.6.0                ReactomePA_1.50.0           cowplot_1.1.3              
 [5] patchwork_1.3.0             BiocManager_1.30.25         apeglm_1.28.0               ggrepel_0.9.6              
 [9] DESeq2_1.46.0               SummarizedExperiment_1.36.0 MatrixGenerics_1.18.0       matrixStats_1.4.1          
[13] GenomicRanges_1.58.0        GenomeInfoDb_1.42.0         RColorBrewer_1.1-3          ggupset_0.4.0              
[17] enrichplot_1.26.3           DOSE_4.0.0                  pheatmap_1.0.12             org.Hs.eg.db_3.20.0        
[21] AnnotationDbi_1.68.0        IRanges_2.40.0              S4Vectors_0.44.0            Biobase_2.66.0             
[25] BiocGenerics_0.52.0         clusterProfiler_4.14.4      lubridate_1.9.4             forcats_1.0.0              
[29] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.5                
[33] tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1               tidyverse_2.0.0  
