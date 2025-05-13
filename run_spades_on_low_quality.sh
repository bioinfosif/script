#!/bin/bash

# Répertoires
TRIMMED_DIR="CHIKV_out/trimmed"
OUTPUT_DIR="CHIKV_out/denovo_assembly"
SAMPLES_FILE="low_quality_samples.txt"

mkdir -p "$OUTPUT_DIR"

echo "🔧 Assemblage de novo avec SPAdes à partir des reads nettoyés..."

while read -r sample; do
    echo "🧬 Assemblage de $sample..."

    R1="${TRIMMED_DIR}/${sample}_R1_trimmed.fastq.gz"
    R2="${TRIMMED_DIR}/${sample}_R2_trimmed.fastq.gz"

    if [[ -f "$R1" && -f "$R2" ]]; then
        spades -t 12 --only-assembler -1 "$R1" -2 "$R2" -o "$OUTPUT_DIR/$sample" --careful
    else
        echo "⚠️ Fichiers trimmed introuvables pour $sample, ignoré."
    fi
done < "$SAMPLES_FILE"

echo "✅ Assemblages de novo terminés dans : $OUTPUT_DIR/"
