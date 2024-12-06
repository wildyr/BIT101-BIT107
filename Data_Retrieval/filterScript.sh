#!/bin/bash

# Input and output file paths
input_file="GSE217504_series_metadata_matrix.txt"
output_file="filtered_metadata.txt"
variables_file="variable_list_edit.txt"

# Ensure the variables list file exists
if [[ ! -f "$variables_file" ]]; then
    echo "Error: $variables_file not found!"
    exit 1
fi

# Read variables of interest into an array
mapfile -t variables_to_keep < "$variables_file"

# Create a pattern for grep
pattern=$(printf "|%s" "${variables_to_keep[@]}")
pattern=${pattern:1}  # Remove the leading '|'

# Extract relevant lines based on the pattern
grep -E "$pattern" "$input_file" > "$output_file"

echo "Metadata filtered from rownames in $variables_file and saved to $output_file"
