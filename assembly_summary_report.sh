#!/bin/bash

# Dossiers
REPORT_DIR="CHIKV_out/assembly_reports"
SUMMARY_FILE="CHIKV_out/assembly_summary.tsv"

# Entête du tableau
echo -e "Sample\tContigs_≥500bp\tTotal_bp\tLargest_contig\tN50" > "$SUMMARY_FILE"

echo "📋 Création du rapport global à partir des rapports individuels..."

# Parcours de tous les rapports individuels
for report in "$REPORT_DIR"/*_report.txt; do
    sample=$(basename "$report" _report.txt)

    # Extraction des valeurs depuis le rapport texte
    contigs=$(grep "Nombre de contigs" "$report" | awk -F: '{gsub(/ /,"",$2); print $2}')
    total_bp=$(grep "Taille totale" "$report" | awk -F: '{gsub(/ /,"",$2); print $2}')
    largest=$(grep "Plus grand contig" "$report" | awk -F: '{gsub(/ bp/,"",$2); gsub(/ /,"",$2); print $2}')
    n50=$(grep "N50" "$report" | awk -F: '{gsub(/ bp/,"",$2); gsub(/ /,"",$2); print $2}')

    # Ajout à la ligne du tableau
    echo -e "${sample}\t${contigs}\t${total_bp}\t${largest}\t${n50}" >> "$SUMMARY_FILE"
done

echo "✅ Rapport global généré : $SUMMARY_FILE"
