#!/bin/bash

ALIGNMENT_DIR="CHIKV_out/alignment"
REPORT="$ALIGNMENT_DIR/alignment_report.csv"

echo "sample,mapped_reads,total_reads,percent_mapped,avg_depth,genome_covered_percent" > "$REPORT"

for flagstat_file in "$ALIGNMENT_DIR"/*_flagstat.txt; do
    sample_name=$(basename "$flagstat_file" _flagstat.txt)
    
    # Lire mapped/total reads
    mapped_line=$(grep 'mapped (' "$flagstat_file" | head -n1)
    mapped_reads=$(echo "$mapped_line" | awk '{print $1}')
    total_reads=$(grep 'in total' "$flagstat_file" | awk '{print $1}')
    percent_mapped=$(echo "$mapped_line" | grep -oP '\(\K[0-9.]+(?=%)')

    # Moyenne de profondeur
    depth_file="$ALIGNMENT_DIR/${sample_name}_depth.txt"
    if [[ -f "$depth_file" ]]; then
        avg_depth=$(awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}' "$depth_file")
    else
        avg_depth="NA"
    fi

    # Couverture
    coverage_file="$ALIGNMENT_DIR/${sample_name}_coverage.txt"
    if [[ -f "$coverage_file" ]]; then
        covered=$(grep -v "^#rname" "$coverage_file" | awk '{print $6}')
    else
        covered="NA"
    fi

    echo "$sample_name,$mapped_reads,$total_reads,$percent_mapped,$avg_depth,$covered" >> "$REPORT"
done

echo "✅ Rapport généré : $REPORT"
