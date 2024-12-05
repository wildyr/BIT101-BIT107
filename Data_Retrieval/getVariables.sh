# Create variable_list.txt file containing all row headers from metadata. We can then manually edit this list for the row headers we want from the metadata file.
grep "^!" GSE217504_series_metadata_matrix.txt | cut -f1 | sort | uniq > variable_list.txt

