#!/bin/bash

# Input and output file paths
input_file="filtered_metadata.txt"
output_file="prepared_metadata.csv"

# Temporary file for intermediate steps
temp_file="temp_metadata.txt"

# Function to clean and rename variable names
clean_variable_name() {
    local name="$1"
    # Remove "!" and replace underscores with spaces
    echo "$name" | sed -e 's/^!//' -e 's/_/ /g' -e 's/ch1/ (ch1)/g'
}

# Prepare the header
{
    # Process each line in the input file
    while IFS=$'\t' read -r first_column rest_of_line; do
        # Clean the variable name
        cleaned_name=$(clean_variable_name "$first_column")
        # Write the cleaned name followed by the rest of the line
        echo -e "$cleaned_name\t$rest_of_line"
    done < "$input_file"
} > "$temp_file"

# Convert the cleaned file to CSV format
# Replace tabs with commas
sed 's/\t/,/g' "$temp_file" > "$output_file"

# Clean up temporary file
rm "$temp_file"

echo "Metadata prepared for R and saved to $output_file"

