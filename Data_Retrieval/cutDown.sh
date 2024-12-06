#!/bin/bash

# Input and output file paths
input_file="prepared_metadata.csv"
output_file="final_metadata.csv"

# Create a new output file by filtering lines that start with "treatment", "hours", or "description"
# and capitalize the first character of each row name.
{
    # Read each line of the input file
    while IFS=, read -r first_col rest_of_line; do
        # Check if the first column starts with "treatment", "hours", or "description"
        if [[ "$first_col" == "treatment" || "$first_col" == "hours" || "$first_col" == "description" ]]; then
            # Capitalize the first character of the first column
            capitalized_first_col=$(echo "$first_col" | sed 's/^\(.\)/\U\1/')
            # Print the capitalized first column with the rest of the line
            echo "$capitalized_first_col,$rest_of_line"
        fi
    done < "$input_file"
} > "$output_file"

echo "Filtered and capitalized rows saved to $output_file"
