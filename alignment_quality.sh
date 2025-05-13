#!/bin/bash

# Répertoire contenant les BAM triés
ALIGNMENT_DIR="CHIKV_out/alignment"

echo "🔍 Évaluation de la qualité des alignements..."

for bam_file in "$ALIGNMENT_DIR"/*_sorted.bam; do
    sample_name=$(basename "$bam_file" _sorted.bam)
    
    echo "📌 Traitement de $sample_name..."

    # Flagstat
    samtools flagstat "$bam_file" > "$ALIGNMENT_DIR/${sample_name}_flagstat.txt"

    # Depth
    samtools depth "$bam_file" > "$ALIGNMENT_DIR/${sample_name}_depth.txt"

    # Coverage (si disponible)
    if samtools coverage --help &> /dev/null; then
        samtools coverage "$bam_file" > "$ALIGNMENT_DIR/${sample_name}_coverage.txt"
    else
        echo "⚠️ samtools coverage non disponible" > "$ALIGNMENT_DIR/${sample_name}_coverage.txt"
    fi
done

echo "✅ Analyse de qualité d’alignement terminée pour tous les échantillons."
