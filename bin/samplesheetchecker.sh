#!/bin/bash

# Function to check if a file exists
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 does not exist."
        exit 1
    fi
}

# Function to validate the samplesheet
validate_samplesheet() {
    samplesheet_path=$1

    # Check if the samplesheet path is correct
    check_file_exists $samplesheet_path

    # Read the header of the samplesheet
    header=$(head -n 1 $samplesheet_path)
    expected_header="sample,fq1,fq2,platform,library,center"

    # Check if the header matches the expected format
    if [ $header != $expected_header ]; then
        echo "Error: Invalid header format."
        exit 1
    fi

    # Check if all rows are CSV and not tab-delimited
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" == *$'\t'* ]]; then
            echo "Error: Found a column with a tab delimiter, expected CSV."
            exit 1
        fi
    done < <(tail -n +2 $samplesheet_path)  # Skip the header line
}

# Main script execution
input_csv="$1"
output_csv="validated_samplesheet.csv"

# Validate samplesheet
validate_samplesheet "$input_csv"

# If validation passes, copy the contents to validated_samplesheet.csv
cat "$input_csv" > "$output_csv"