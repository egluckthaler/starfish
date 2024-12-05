#!/bin/bash

# Utility script for combining filt and filt_intersect output from starfish annotate
# across multiple independent runs

# Function to display help menu
display_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Utility script for combining filt and filt_intersect output from starfish annotate"
    echo ""
    echo "Options:"
    echo "  -i, --input FILE        Input TSV file with genome codes in the first field (required)"
    echo "  -a, --analysis PREFIX   Path to root analysis directory containing an output directory for each genome (required)"
    echo "  -o, --output PREFIX     Path and Prefix for output files (required)"
    echo "  -h, --help              Display this help menu"
    echo ""
    echo "Example:"
    echo "  $0 -i ome2assembly.txt -a starfish_run1 -o all_annotations"
}

# Initialize variables
INPUT_FILE=""
ANALYSIS_PREFIX=""
OUTPUT_PREFIX=""
HELP=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -a|--analysis)
            ANALYSIS_PREFIX="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -h|--help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            display_help
            exit 1
            ;;
    esac
done

# Display help or validate required arguments
if [[ "$HELP" = true ]] || [ -z "$INPUT_FILE" ] || [ -z "$ANALYSIS_PREFIX" ] || [ -z "$OUTPUT_PREFIX" ]; then
    display_help
    exit 1
fi

output_gff=${OUTPUT_PREFIX}.gff
output_ids=${OUTPUT_PREFIX}.ids
output_fas=${OUTPUT_PREFIX}.fas

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

# Initialize counters
count_filt_intersect=0
count_filt=0
count_none=0

# Clear the output file or create it if it doesn't exist
> "$output_gff"
> "$output_ids"
> "$output_fas"

# Read each genome code from the input file
while IFS=$'\t' read -r genome_code rest_of_line || [ -n "$genome_code" ]; do
    
    # Trim any potential whitespace from the genome code
    genome_code=$(echo "$genome_code" | xargs)

    # Define the directory name based on the genome code
    dir_name="${ANALYSIS_PREFIX}/$genome_code"

    # Check for 'filt_intersect' file first
    if [ -f "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt_intersect.ids" ]; then
        cat "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt_intersect.ids" >> "$output_ids"
        cat "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt_intersect.fas" >> "$output_fas"
        cat "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt_intersect.gff" >> "$output_gff"
        ((count_filt_intersect++))
    elif [ -f "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt.ids" ]; then
        # If 'filt_intersect' doesn't exist, use 'filt' file
        cat "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt.ids" >> "$output_ids"
        cat "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt.fas" >> "$output_fas"
        cat "${dir_name}/${genome_code}.${ANALYSIS_PREFIX}.filt.gff" >> "$output_gff"
        ((count_filt++))
    else
        echo "No 'filt_intersect' or 'filt' file found for genome code: $genome_code"
        ((count_none++))
    fi
done < "$INPUT_FILE"

echo "Total genomes with 'filt_intersect' files: $count_filt_intersect"
echo "Total genomes with 'filt' files: $count_filt"
echo "Total genomes with no files: $count_none"
echo "Processing complete"