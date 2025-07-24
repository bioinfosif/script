#!/bin/bash

# Prokka batch processing script
# This script runs Prokka annotation on all FASTA files in the genomes/ directory

# Set up directories
GENOME_DIR="genomes"
OUTPUT_DIR="prokka_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if Prokka is installed
if ! command -v prokka &> /dev/null; then
    echo "Error: Prokka is not installed or not in PATH"
    echo "Please install Prokka first: conda install -c bioconda prokka"
    exit 1
fi

# Function to get base name without extension
get_basename() {
    local filename="$1"
    echo "${filename%.*}"
}

# Counter for progress tracking
total_files=$(ls -1 "$GENOME_DIR"/*.fasta 2>/dev/null | wc -l)
current=0

echo "Found $total_files FASTA files to process"
echo "Starting Prokka annotation..."

# Process each FASTA file
for fasta_file in "$GENOME_DIR"/*.fasta; do
    if [ -f "$fasta_file" ]; then
        current=$((current + 1))
        
        # Get the base name for output directory
        filename=$(basename "$fasta_file")
        base_name=$(get_basename "$filename")
        
        echo "[$current/$total_files] Processing: $filename"
        
        # Create output subdirectory for this genome
        output_subdir="$OUTPUT_DIR/$base_name"
        
        # Run Prokka
        prokka \
            --outdir "$output_subdir" \
            --prefix "$base_name" \
            --kingdom Bacteria \
            --locustag "$base_name" \
            --cpus 30 \
            --force \
            "$fasta_file"
        
        if [ $? -eq 0 ]; then
            echo "✓ Successfully annotated $filename"
        else
            echo "✗ Failed to annotate $filename"
        fi
        
        echo "---"
    fi
done

echo "Prokka annotation completed!"
echo "Results are in: $OUTPUT_DIR/"
echo ""
echo "Output files for each genome include:"
echo "  *.gff - Annotation in GFF3 format"
echo "  *.gbk - Genbank format"
echo "  *.fna - Nucleotide FASTA file"
echo "  *.faa - Protein FASTA file"
echo "  *.ffn - Feature nucleotide sequences"
echo "  *.sqn - Sequin file for GenBank submission"
echo "  *.fsa - Nucleotide FASTA file for GenBank submission"
echo "  *.tbl - Feature table for GenBank submission"
echo "  *.err - Unacceptable annotations"
echo "  *.log - Log file"
echo "  *.txt - Statistics"
