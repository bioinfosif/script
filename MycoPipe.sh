#!/bin/bash

# Pipeline d'analyse bioinformatique pour Mycobacterium
# Basé sur la méthodologie décrite dans l'article

# Configuration des variables
OUTPUT_DIR="mycobacterium_analysis"
THREADS=16
MIN_SCAFFOLD_LENGTH=800
DEPTH_THRESHOLD=0.25  # 25% de la profondeur moyenne

# Créer la structure de répertoires
mkdir -p $OUTPUT_DIR/{1_fastqc,2_trimmed,3_assembly,4_quality,5_taxonomy,6_resistance,7_phylogeny}

# Fonction pour traiter un échantillon
process_sample() {
    sample_name=$1
    r1=$2
    r2=$3
    
    echo "Traitement de l'échantillon: $sample_name"
    
    # 1. Contrôle qualité avec FastQC
    echo "Étape 1: Contrôle qualité avec FastQC"
    fastqc -o $OUTPUT_DIR/1_fastqc -t $THREADS $r1 $r2
    
    # 2. Trimming avec Trimmomatic
    echo "Étape 2: Trimming avec Trimmomatic"
    trimmomatic PE -threads $THREADS $r1 $r2 \
        $OUTPUT_DIR/2_trimmed/${sample_name}_R1_paired.fastq.gz $OUTPUT_DIR/2_trimmed/${sample_name}_R1_unpaired.fastq.gz \
        $OUTPUT_DIR/2_trimmed/${sample_name}_R2_paired.fastq.gz $OUTPUT_DIR/2_trimmed/${sample_name}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # 3. Assemblage avec Unicycler
    echo "Étape 3: Assemblage avec Unicycler"
    unicycler -1 $OUTPUT_DIR/2_trimmed/${sample_name}_R1_paired.fastq.gz \
              -2 $OUTPUT_DIR/2_trimmed/${sample_name}_R2_paired.fastq.gz \
              -o $OUTPUT_DIR/3_assembly/$sample_name \
              -t $THREADS
    
    # Filtrage des scaffolds courts et des contaminants potentiels
    echo "Filtrage des scaffolds"
    python3 filter_scaffolds.py \
        --input $OUTPUT_DIR/3_assembly/$sample_name/assembly.fasta \
        --output $OUTPUT_DIR/3_assembly/$sample_name/filtered_assembly.fasta \
        --min_length $MIN_SCAFFOLD_LENGTH \
        --depth_threshold $DEPTH_THRESHOLD
    
    # 4. Évaluation de la qualité du génome avec CheckM
    echo "Étape 4: Évaluation de la qualité avec CheckM"
    checkm lineage_wf -t $THREADS -x fasta $OUTPUT_DIR/3_assembly/$sample_name/ \
        $OUTPUT_DIR/4_quality/$sample_name
    
    # 5. Analyse de la taxonomie avec TYGS (nécessite une soumission manuelle)
    echo "Étape 5: Préparation des fichiers pour TYGS"
    cp $OUTPUT_DIR/3_assembly/$sample_name/filtered_assembly.fasta $OUTPUT_DIR/5_taxonomy/${sample_name}.fasta
    echo "Soumettre $OUTPUT_DIR/5_taxonomy/${sample_name}.fasta au serveur TYGS: https://tygs.dsmz.de/"
    
    # 6. Analyse des résistances et détermination des lignées avec TB-Profiler
    echo "Étape 6: Analyse des résistances avec TB-Profiler"
    tb-profiler profile -1 $OUTPUT_DIR/2_trimmed/${sample_name}_R1_paired.fastq.gz \
                       -2 $OUTPUT_DIR/2_trimmed/${sample_name}_R2_paired.fastq.gz \
                       -p $sample_name \
                       -t $THREADS \
                       --dir $OUTPUT_DIR/6_resistance
}

# Fonction pour l'analyse phylogénétique de tous les échantillons
run_phylogenetic_analysis() {
    echo "Étape 7: Alignement des séquences avec Mugsy"
    # Créer un fichier de liste des génomes
    find $OUTPUT_DIR/3_assembly -name "filtered_assembly.fasta" > $OUTPUT_DIR/7_phylogeny/genome_list.txt
    
    # Exécuter Mugsy pour l'alignement
    mugsy --directory $OUTPUT_DIR/7_phylogeny \
          --prefix mycobacterium_alignment \
          $(cat $OUTPUT_DIR/7_phylogeny/genome_list.txt)
    
    # Convertir l'alignement Mugsy en format FASTA pour IQ-TREE
    python3 convert_mugsy_to_fasta.py \
        --input $OUTPUT_DIR/7_phylogeny/mycobacterium_alignment.maf \
        --output $OUTPUT_DIR/7_phylogeny/mycobacterium_alignment.fasta
    
    echo "Construction de l'arbre phylogénétique avec IQ-TREE2"
    # Exécuter IQ-TREE2 avec le modèle GTR et 1000 bootstraps
    iqtree2 -s $OUTPUT_DIR/7_phylogeny/mycobacterium_alignment.fasta \
           -m GTR \
           -b 1000 \
           -T $THREADS
    
    echo "Analyse phylogénétique terminée. L'arbre peut être visualisé avec iTOL."
}

# Rendre les scripts Python exécutables
chmod +x filter_scaffolds.py convert_mugsy_to_fasta.py

# Usage du pipeline
echo "=== PIPELINE D'ANALYSE GÉNOMIQUE POUR MYCOBACTERIUM ==="
echo "Pour traiter un échantillon individuel:"
echo "  process_sample nom_echantillon read1.fastq.gz read2.fastq.gz"
echo ""
echo "Pour lancer l'analyse phylogénétique après traitement des échantillons:"
echo "  run_phylogenetic_analysis"
echo ""
echo "Exemple de pipeline complet:"
echo '  # Traiter chaque échantillon'
echo '  process_sample sample1 data/sample1_R1.fastq.gz data/sample1_R2.fastq.gz'
echo '  process_sample sample2 data/sample2_R1.fastq.gz data/sample2_R2.fastq.gz'
echo "  # Lancer l'analyse phylogénétique"
echo '  run_phylogenetic_analysis'
