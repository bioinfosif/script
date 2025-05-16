#!/bin/bash

# --- CONFIGURATION ---
READ_DIR="data"
OUT_DIR="results_tbprofiler"
THREADS=28
ADAPTERS="/home/genomic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
TRIMMOMATIC_JAR="~/Trimmomatic-0.39/trimmomatic-0.39.jar"

mkdir -p "$OUT_DIR"

# --- BOUCLE SUR LES ÉCHANTILLONS ---
for R1 in ${READ_DIR}/*_R1_001.fastq.gz; do
    # Récupère le nom de base
    BASENAME=$(basename "$R1")
    
    # Remplace _R1_001 par _R2_001 pour retrouver le bon fichier pairé
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

    # Récupère le nom de l’échantillon sans l’extension
    SAMPLE=$(echo "$BASENAME" | cut -d"_" -f1)

    echo "[ $SAMPLE ] Traitement en cours..."

    # --- 1. Trimming avec Trimmomatic ---
    TRIM_DIR="$OUT_DIR/${SAMPLE}/trimmed"
    mkdir -p "$TRIM_DIR"
    TRIM_R1="${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    TRIM_R2="${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    UNPAIRED1="${TRIM_DIR}/unpaired1.fq.gz"
    UNPAIRED2="${TRIM_DIR}/unpaired2.fq.gz"

    java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $THREADS \
        "$R1" "$R2" \
        "$TRIM_R1" "$UNPAIRED1" \
        "$TRIM_R2" "$UNPAIRED2" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        SLIDINGWINDOW:4:20 MINLEN:50

    # --- 2. TB-Profiler ---
    TB_OUT="$OUT_DIR/${SAMPLE}"
    mkdir -p "$TB_OUT"
    cd "$TB_OUT"
    pwd
    tb-profiler profile \
        --read1 "../../$TRIM_R1" \
        --read2 "../../$TRIM_R2" \
        --threads $THREADS \
        --prefix "$SAMPLE" \
        --csv

    cd - > /dev/null
    pwd
done

echo "[✔] Analyse terminée. Résultats dans '$OUT_DIR'"
