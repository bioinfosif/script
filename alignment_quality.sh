#!/bin/bash

# BAM files Directory
ALIGNMENT_DIR="CHIKV_out/alignment"

echo "Alignment quality processing..."

for bam_file in "$ALIGNMENT_DIR"/*_sorted.bam; do
    sample_name=$(basename "$bam_file" _sorted.bam)
    
    echo "Processing $sample_name..."

    # Flagstat
    samtools flagstat "$bam_file" > "$ALIGNMENT_DIR/${sample_name}_flagstat.txt"

    # Depth
    samtools depth "$bam_file" > "$ALIGNMENT_DIR/${sample_name}_depth.txt"

    # Coverage (if possible)
    if samtools coverage --help &> /dev/null; then
        samtools coverage "$bam_file" > "$ALIGNMENT_DIR/${sample_name}_coverage.txt"
    else
        echo "samtools coverage not available" > "$ALIGNMENT_DIR/${sample_name}_coverage.txt"
    fi
done
 
echo "Alignment quality analysis completed for all samples."
