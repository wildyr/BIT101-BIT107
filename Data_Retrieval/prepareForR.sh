#!/bin/bash

# File paths
input_file="converted_metadata.csv"
output_file="prepared_metadata.csv"

# Define the new row headers
new_headers=("cell line" "cell type" "genotype" "treatment" "hours" "description")

# Read the input file line by line and replace the row header
{
    # Initialize a counter for the new headers
    i=0

    # Read each line of the input file
    while IFS=, read -r first_col rest_of_line; do
        # Replace the first column with the new header if we have one
        if [ $i -lt ${#new_headers[@]} ]; then
            echo "${new_headers[$i]},$rest_of_line"
            ((i++))
        else
            # If we run out of new headers, just print the line as it is
            echo "$first_col,$rest_of_line"
        fi
    done < "$input_file"
} > "$output_file"

# Loop through each header and remove instances where it appears followed by a colon
for header in "${new_headers[@]}"; do
    # Use sed to remove the header followed by a colon from the entire document
    sed -i "s/\b$header\b: //g" "$output_file"
done

sed -i 's/time point_(in_hours): //g' "$output_file"

echo "Row headers updated and saved to $output_file"

