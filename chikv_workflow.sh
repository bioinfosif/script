
#!/bin/bash

# Script pour le traitement complet des échantillons CHIKV d'intérêt
# Auteur: Claude
# Date: 21 mai 2025

# ===========================
# Paramètres à configurer
# ===========================

# Liste des échantillons d'intérêt (>15% CHIKV)
# Échantillons avec >15% de CHIKV identifiés
SAMPLES=(
    "13_S6" "45_S28" "47_S29" "27_S15" "19_S9" 
    "3_S2" "39_S24" "51_S31" "9_S3" "20_S10" 
    "50_S30" "2_S1" "38_S23" "34_S20" "31_S18" 
    "40_S25" "37_S22"
)

# Chemin vers les fichiers FASTQ bruts
RAW_DATA_DIR="../data"

# Nombre de threads à utiliser
THREADS=10

# Référence CHIKV (URL NCBI)
REF_URL="https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=NC_004162.2"

# Paramètres Trimmomatic
TRIM_PARAMS="ILLUMINACLIP:/home/biobacteria/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50"

# Seuils de filtrage
CONTAMINATION_THRESHOLD=5.0
COMPLETENESS_THRESHOLD=90.0
HETEROGENEITY_THRESHOLD=10.0

# ===========================
# Création des répertoires
# ===========================

echo "Création des répertoires de travail..."
mkdir -p fastqc_results trimmed_reads references alignments consensus quast_results checkm2_results checkm2_input

# ===========================
# 1. Contrôle qualité et trimming
# ===========================

echo "====== ÉTAPE 1: Contrôle qualité et trimming ======"

# 1.1 Analyse FastQC
echo "Exécution de FastQC sur les échantillons bruts..."
for sample in "${SAMPLES[@]}"; do
    echo "  Processing ${sample}..."
    fastqc -t $THREADS -o fastqc_results ${RAW_DATA_DIR}/${sample}_L001_R1_001.fastq.gz ${RAW_DATA_DIR}/${sample}_L001_R2_001.fastq.gz
done

# Rapport agrégé avec MultiQC
echo "Génération du rapport MultiQC..."
multiqc fastqc_results -o fastqc_results

# 1.2 Trimming avec Trimmomatic
echo "Exécution de Trimmomatic..."
for sample in "${SAMPLES[@]}"; do
    echo "  Trimming ${sample}..."
    R1=${RAW_DATA_DIR}/${sample}_L001_R1_001.fastq.gz
    R2=${RAW_DATA_DIR}/${sample}_L001_R2_001.fastq.gz
    
    java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $THREADS \
        $R1 $R2 \
        trimmed_reads/${sample}_R1_paired.fastq.gz trimmed_reads/${sample}_R1_unpaired.fastq.gz \
        trimmed_reads/${sample}_R2_paired.fastq.gz trimmed_reads/${sample}_R2_unpaired.fastq.gz \
        $TRIM_PARAMS
done

# Vérification post-trimming
echo "Vérification de la qualité après trimming..."
mkdir -p fastqc_trimmed
for sample in "${SAMPLES[@]}"; do
    fastqc -t $THREADS -o fastqc_trimmed trimmed_reads/${sample}_R1_paired.fastq.gz trimmed_reads/${sample}_R2_paired.fastq.gz
done

multiqc fastqc_trimmed -o fastqc_trimmed

# ===========================
# 2. Assemblage basé sur référence
# ===========================

echo "====== ÉTAPE 2: Assemblage basé sur référence ======"


# 2.1 Téléchargement de la référence
echo "Téléchargement de la référence CHIKV..."
wget -O references/CHIKV_ref.fasta "$REF_URL"


# 2.2 Indexation de la référence avec BWA
echo "Indexation de la référence..."
bwa index references/CHIKV_ref.fasta
samtools faidx references/CHIKV_ref.fasta

# 2.3 Alignement des lectures
echo "Alignement des lectures sur la référence..."
for sample in "${SAMPLES[@]}"; do
    echo "  Alignement de ${sample}..."
    
    # Alignement BWA
    bwa mem -t $THREADS references/CHIKV_ref.fasta \
        trimmed_reads/${sample}_R1_paired.fastq.gz \
        trimmed_reads/${sample}_R2_paired.fastq.gz \
        > alignments/${sample}.sam
    
    # Conversion et traitement du BAM
    echo "  Traitement BAM pour ${sample}..."
    samtools view -bS alignments/${sample}.sam > alignments/${sample}.bam
    samtools sort alignments/${sample}.bam -o alignments/${sample}.sorted.bam
    samtools index alignments/${sample}.sorted.bam
    
    # Nettoyage
    rm alignments/${sample}.sam alignments/${sample}.bam
    
    # Statistiques d'alignement
    samtools flagstat alignments/${sample}.sorted.bam > alignments/${sample}.flagstat
    samtools depth -a alignments/${sample}.sorted.bam > alignments/${sample}.depth
done

# 2.4 Génération des consensus
echo "Génération des consensus..."
for sample in "${SAMPLES[@]}"; do
    echo "  Génération consensus pour ${sample}..."
    
    # Appel de variants
    bcftools mpileup -f references/CHIKV_ref.fasta alignments/${sample}.sorted.bam | \
    bcftools call -mv -Oz -o consensus/${sample}.vcf.gz
    
    # Indexation du VCF
    bcftools index consensus/${sample}.vcf.gz
    
    # Génération du consensus
    bcftools consensus -f references/CHIKV_ref.fasta consensus/${sample}.vcf.gz > consensus/${sample}.consensus.fasta
    
    # Renommage du consensus
    sed -i "s/>.*/>$sample/" consensus/${sample}.consensus.fasta
    
    # Copier le consensus pour l'analyse CheckM2
    cp consensus/${sample}.consensus.fasta checkm2_input/
done

# ===========================
# 3. Évaluation de qualité
# ===========================

echo "====== ÉTAPE 3: Évaluation de qualité ======"

# 3.1 Analyse QUAST
echo "Exécution de QUAST..."
quast consensus/*.consensus.fasta \
    -r references/CHIKV_ref.fasta \
    -o quast_results \
    --threads $THREADS

# 3.2 Analyse CheckM2
echo "Exécution de CheckM2..."
checkm2 predict --threads $THREADS \
    --input checkm2_input/* \
    --output-directory checkm2_results \
    --force

# ===========================
# 4. Filtrage et sélection
# ===========================

echo "====== ÉTAPE 4: Filtrage et sélection des meilleurs consensus ======"

# Exécuter le script de filtrage
echo "Filtrage des consensus..."
export CONTAMINATION_THRESHOLD=$CONTAMINATION_THRESHOLD
export COMPLETENESS_THRESHOLD=$COMPLETENESS_THRESHOLD
export HETEROGENEITY_THRESHOLD=$HETEROGENEITY_THRESHOLD

python3 filter_consensus.py


# ===========================
# 5. Visualisation (optionnel)
# ===========================

echo "====== ÉTAPE 5: Génération des visualisations ======"

# Génération de graphiques de couverture
echo "Génération des graphiques de couverture..."
for sample in "${SAMPLES[@]}"; do
    echo "  Graphique pour ${sample}..."
    
    # Créer un script R pour la visualisation
    cat > plot_coverage_${sample}.R << EOL
library(ggplot2)

# Lire les données de couverture
coverage_data <- read.table("alignments/${sample}.depth", 
                            col.names=c("Chromosome", "Position", "Depth"))

# Créer le graphique
p <- ggplot(coverage_data, aes(x=Position, y=Depth)) +
    geom_line(color="blue") +
    labs(title="Profondeur de couverture pour ${sample}",
         x="Position sur le génome", 
         y="Profondeur") +
    theme_minimal()

# Enregistrer le graphique
ggsave("alignments/${sample}_coverage.png", plot=p, width=10, height=6, dpi=300)
EOL

    # Exécuter le script R si R est installé
    if command -v Rscript &> /dev/null; then
        Rscript plot_coverage_${sample}.R
    else
        echo "R n'est pas installé. Graphique de couverture non généré."
    fi
done

# ===========================
# Résumé final
# ===========================

echo "====== RÉSUMÉ FINAL ======"
echo "Traitement terminé!"
echo "Résultats disponibles dans:"
echo "  - Rapports FastQC: ./fastqc_results"
echo "  - Lectures filtrées: ./trimmed_reads"
echo "  - Alignements: ./alignments"
echo "  - Consensus: ./consensus"
echo "  - Rapport QUAST: ./quast_results"
echo "  - Rapport CheckM2: ./checkm2_results"
echo "  - Meilleurs échantillons: ./filtered_samples.csv"
echo "  - Meilleur consensus: ./best_consensus.fasta"
echo ""
echo "Pour visualiser le meilleur consensus:"
echo "less best_consensus.fasta"
echo ""
echo "Pour visualiser les statistiques de filtrage:"
echo "cat filtered_samples.csv"
