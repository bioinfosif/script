#!/bin/bash
# Script complet pour méta-assemblage avec Flye

# Variables
READS="data/AgroC1.fastq.gz"
OUTPUT_DIR="AgroC1"
THREADS=10
GENOME_SIZE="100m"

# Étape 1: Vérifier les données d'entrée
echo "=== Vérification des données ==="
if [[ ! -f "$READS" ]]; then
    echo "Erreur: Fichier de lectures $READS introuvable"
    exit 1
fi

echo "Nombre de lectures: $(wc -l < $READS | awk '{print $1/4}')"
echo "Taille du fichier: $(ls -lh $READS | awk '{print $5}')"

# Étape 2: Lancer Flye
echo "=== Lancement de Flye ==="
flye --meta \
     --nano-raw "$READS" \
     --out-dir "$OUTPUT_DIR" \
     --threads "$THREADS" \
     

# Étape 3: Vérifier les résultats
echo "=== Résultats ==="
if [[ -f "$OUTPUT_DIR/assembly.fasta" ]]; then
    echo "Assemblage réussi!"
    echo "Fichier de sortie: $OUTPUT_DIR/assembly.fasta"
    echo "Statistiques:"
    grep -c ">" "$OUTPUT_DIR/assembly.fasta"
    echo "contigs générés"
else
    echo "Erreur: Assemblage échoué"
    exit 1
fi

# =============================================================================
# ANALYSE POST-ASSEMBLAGE
# =============================================================================

# Statistiques de l'assemblage
echo "=== Statistiques de l'assemblage ==="

# Nombre de contigs
CONTIGS=$(grep -c ">" "$OUTPUT_DIR/assembly.fasta")
echo "Nombre de contigs: $CONTIGS"

# Taille totale de l'assemblage
TOTAL_SIZE=$(awk '/^>/ {next} {total += length($0)} END {print total}' "$OUTPUT_DIR/assembly.fasta")
echo "Taille totale: $TOTAL_SIZE bp"

# N50 (nécessite un script séparé ou un outil comme assembly-stats)
# assembly-stats $OUTPUT_DIR/assembly.fasta

# =============================================================================
# COMMANDES DE MONITORING
# =============================================================================

# Surveiller l'avancement (dans un autre terminal)
# tail -f $OUTPUT_DIR/flye.log

# Vérifier l'utilisation des ressources
# htop
# ou
# watch -n 5 'ps aux | grep flye'

# =============================================================================
# GESTION DES ERREURS COURANTES
# =============================================================================

# Si erreur de mémoire insuffisante:
# 1. Réduire --genome-size
# 2. Réduire --threads
# 3. Utiliser --asm-coverage pour limiter la couverture

# Si assemblage trop fragmenté:
# 1. Augmenter --min-overlap
# 2. Augmenter --iterations
# 3. Vérifier la qualité des données d'entrée

# Si temps d'exécution trop long:
# 1. Augmenter --threads (si possible)
# 2. Réduire --iterations
# 3. Utiliser --genome-size plus précis
