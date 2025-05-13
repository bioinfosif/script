#!/bin/bash

# Fichiers d'entrée
FASTA="gisaid_sequences.fasta"
IDS="ID_selected.txt"
OUTPUT="CHIKV_selected.fasta"

# Vérifier que seqkit est installé
if ! command -v seqkit &> /dev/null; then
    echo "❌ Erreur : 'seqkit' n'est pas installé. Installe-le avec :"
    echo "   conda install -c bioconda seqkit"
    echo "   ou"
    echo "   brew install seqkit"
    exit 1
fi

# Exécuter l'extraction
seqkit grep -f "$IDS" "$FASTA" > "$OUTPUT"

# Vérification
COUNT=$(grep -c "^>" "$OUTPUT")
echo "✅ Extraction terminée : $COUNT séquences extraites dans '$OUTPUT'"
