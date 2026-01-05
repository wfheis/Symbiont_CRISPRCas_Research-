#!/bin/bash

INPUT_FILE="csy4_homolog_IDs_large.txt"
OUTPUT_FILE="csy4F_prot_large.fasta"

# Ensure input file exists
[[ ! -f "$INPUT_FILE" ]] && { echo "Error: $INPUT_FILE not found!"; exit 1; }

# Convert file to Unix format (remove hidden Windows line endings)
sed -i 's/\r$//' "$INPUT_FILE"

# Ensure proper `OR` formatting
QUERY=$(awk '{printf "%s OR ", $1}' "$INPUT_FILE")
QUERY=${QUERY% OR }  # Remove the trailing ' OR'

# Debugging: Print query to verify formatting
echo "Batch Query: $QUERY"

# If the query is empty, exit
if [[ -z "$QUERY" ]]; then
    echo "Error: Query is empty! Check input file formatting."
    exit 1
fi

# Perform batch search and fetch sequences
echo "Fetching sequences for batch query..."
esearch -db protein -query "$QUERY" | efetch -format fasta > "$OUTPUT_FILE"

# Verify the output file
if [[ -s "$OUTPUT_FILE" ]]; then
    echo "✅ All protein sequences saved in $OUTPUT_FILE"
else
    echo "⚠️ No sequences retrieved! Check query formatting or NCBI availability."
fi
