#!/bin/bash

# Répertoires
ASSEMBLY_DIR="CHIKV_out/denovo_assembly"
REPORT_DIR="CHIKV_out/assembly_reports"
QUAST_DIR="CHIKV_out/quast_reports"
MIN_LENGTH=500  # Longueur minimale de contig à conserver

mkdir -p "$REPORT_DIR" "$QUAST_DIR"

echo "📊 Génération des rapports d'assemblage + QUAST + filtrage des contigs < $MIN_LENGTH bp..."

for sample_dir in "$ASSEMBLY_DIR"/*; do
    sample=$(basename "$sample_dir")
    contigs_file="$sample_dir/contigs.fasta"
    filtered_contigs="$sample_dir/contigs_filtered_${MIN_LENGTH}bp.fasta"
    report_file="$REPORT_DIR/${sample}_report.txt"
    quast_sample_dir="$QUAST_DIR/$sample"

    if [[ -f "$contigs_file" ]]; then
        echo "🧪 Traitement de $sample..."

        # Filtrage des contigs >= MIN_LENGTH
        awk -v min=$MIN_LENGTH '
            BEGIN {RS=">"; ORS=""}
            NR > 1 {
                header = $1
                seq = ""
                for (i = 2; i <= NF; i++) seq = seq $i
                gsub("\n", "", seq)
                if (length(seq) >= min) {
                    print ">" header "\n" seq "\n"
                }
            }
        ' "$contigs_file" > "$filtered_contigs"

        # Statistiques manuelles
        total_contigs=$(grep -c "^>" "$filtered_contigs")
        total_bases=$(grep -v "^>" "$filtered_contigs" | tr -d '\n' | wc -c)
        largest_contig=$(awk '/^>/ {if (seq) print length(seq); seq=""} !/^>/ {seq = seq $0} END {print length(seq)}' "$filtered_contigs" | sort -nr | head -n1)
        n50=$(awk '/^>/ {if (seq) {print length(seq)}; seq=""} !/^>/ {seq=seq $0} END {if (seq) print length(seq)}' "$filtered_contigs" | \
            sort -nr | awk -v sum=0 -v total=0 '
                {L[NR]=$1; total+=L[NR]}
                END {
                    half=total/2
                    for(i=1;i<=NR;i++){
                        sum+=L[i]
                        if(sum >= half){
                            print L[i]
                            break
                        }
                    }
                }')

        # Rapport texte rapide
        {
            echo "🧬 Rapport d'assemblage pour : $sample"
            echo "Nombre de contigs >= ${MIN_LENGTH} bp : $total_contigs"
            echo "Taille totale (bp) : $total_bases"
            echo "Plus grand contig : $largest_contig bp"
            echo "N50 : $n50 bp"
        } > "$report_file"



    else
        echo "⚠️ Contigs introuvables pour $sample, ignoré."
    fi
done


