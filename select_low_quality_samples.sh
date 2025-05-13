#!/bin/bash

REPORT="CHIKV_out/alignment/alignment_report.csv"
OUTPUT="low_quality_samples.txt"

# Seuils de qualité
COVERAGE_THRESHOLD=80
MAPPED_THRESHOLD=60
DEPTH_THRESHOLD=600


echo "🔍 Sélection des échantillons de mauvaise qualité..." > "$OUTPUT"
echo "# Critères : mapped<$MAPPED_THRESHOLD%, coverage<$COVERAGE_THRESHOLD%, depth<$DEPTH_THRESHOLD" >> "$OUTPUT"

tail -n +2 "$REPORT" | while IFS=',' read -r sample mapped total percent_mapped avg_depth genome_covered; do
    pm=${percent_mapped%%.*}
    cov=${genome_covered%%.*}
    depth=${avg_depth%%.*}

    if [[ "$pm" -lt $MAPPED_THRESHOLD || "$cov" -lt $COVERAGE_THRESHOLD || "$depth" -lt $DEPTH_THRESHOLD ]]; then
        echo "$sample" >> "$OUTPUT"
    fi
done

echo "✅ Liste des échantillons à assembler (de novo) : $OUTPUT"

