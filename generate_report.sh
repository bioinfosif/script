#!/bin/bash

OUTPUT_DIR="CHIKV_out"
REPORT="$OUTPUT_DIR/final_report.csv"

echo "Sample,Consensus_Length,Total_Variants,N50,Coverage" > "$REPORT"

for fasta in "$OUTPUT_DIR"/consensus/*.fasta; do
    sample=$(basename "$fasta" _consensus.fasta)

    consensus_length=$(grep -v ">" "$fasta" | tr -d '\n' | wc -c)

    vcf_gz="$OUTPUT_DIR/variants/${sample}_variants.vcf.gz"
    if [[ -f "$vcf_gz" ]]; then
        variant_count=$(bcftools view "$vcf_gz" | grep -vc "^#")
    else
        variant_count="NA"
    fi

    quast_report="$OUTPUT_DIR/quast/${sample}/report.txt"
    if [[ -f "$quast_report" ]]; then
        N50=$(grep -P "^N50" "$quast_report" | awk '{print $2}')
        coverage=$(grep -P "^Genome fraction" "$quast_report" | awk '{print $3}')
    else
        N50="NA"
        coverage="NA"
    fi

    echo "$sample,$consensus_length,$variant_count,$N50,$coverage" >> "$REPORT"
done

echo "Rapport généré : $REPORT"

