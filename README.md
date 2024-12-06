# assignment2
Assignment for BIT101 and BIT107 - Short communication style mini-project report

Data Retrieval
A host of scripts to filter the metadata file GSE217504_series_metadata_matrix.txt to the relevant data for analysis.

getVariables.sh - Retrieves all rownames and prints them to variable_list.txt
filterScript.sh - Looks in variable_list_edit.txt and filters our data to those specified - filtered_metadata.txt
NOTE: variable_list_edit.txt has been manually edited for this dataset. If running all scripts consecutively this file and the original metadata file are required.
convert.CSV.sh - Cleans row names up for better readability and converts the file to a csv - converted_metadata.csv
prepareForR.sh - Rewrites duplicate row names and removes superflous information from our data points - prepared_metadata.csv
cutDown.sh - Further filters applied to remove unecessary fields before importing to R for analysis - final_metadata.csv

For ease of use I've also added a pipline.sh to run this scripts consequtively. Please note that both the original metadata txt and the variable_list_edit.txt are required.

These scripts were written for the purpose of analysing this specific dataset, and would require significant changes for different datasets.

