#!/usr/bin/env python3

import os
import pandas as pd
import shutil
import glob

# Charger les résultats de QUAST
quast_report = "quast_results/report.tsv"
quast_df = pd.read_csv(quast_report, sep='\t')

# Charger les résultats de CheckM2
checkm2_report = "checkm2_results/quality_report.tsv"
checkm2_df = pd.read_csv(checkm2_report, sep='\t')

# Récupérer les seuils depuis les variables d'environnement
contamination_threshold = float(os.environ.get('CONTAMINATION_THRESHOLD', 5.0))
completeness_threshold = float(os.environ.get('COMPLETENESS_THRESHOLD', 90.0))
heterogeneity_threshold = float(os.environ.get('HETEROGENEITY_THRESHOLD', 10.0))

print(f"Filtrage avec les seuils suivants:")
print(f"  - Contamination: ≤ {contamination_threshold}%")
print(f"  - Complétude: ≥ {completeness_threshold}%")
print(f"  - Hétérogénéité: ≤ {heterogeneity_threshold} mismatches/100kbp")

# Filtrer les résultats
filtered_samples = []
sample_files = glob.glob("consensus/*.consensus.fasta")

for sample_file in sample_files:
    sample_name = os.path.basename(sample_file).replace('.consensus.fasta', '')
    
    # Chercher dans CheckM2
    checkm2_row = checkm2_df[checkm2_df['Name'].str.contains(sample_name)]
    if checkm2_row.empty:
        print(f"Attention: {sample_name} non trouvé dans les résultats CheckM2")
        continue
        
    contamination = checkm2_row['Contamination'].iloc[0]
    completeness = checkm2_row['Completeness'].iloc[0]
    
    # Chercher dans QUAST
    quast_row = quast_df[quast_df['Assembly'].str.contains(sample_name)]
    if quast_row.empty:
        print(f"Attention: {sample_name} non trouvé dans les résultats QUAST")
        continue
        
    n50 = quast_row['N50'].iloc[0]
    genome_fraction = quast_row['Genome fraction (%)'].iloc[0]
    
    # Vérifier si les colonnes existent
    if 'mismatches per 100 kbp' in quast_row.columns:
        mismatches = quast_row['mismatches per 100 kbp'].iloc[0]
    else:
        # Valeur par défaut si la colonne n'existe pas
        mismatches = float('inf')
        print(f"Attention: Colonne 'mismatches per 100 kbp' non trouvée pour {sample_name}")
    
    # Vérifier si l'échantillon répond aux critères
    if (contamination <= contamination_threshold and
        completeness >= completeness_threshold and
        mismatches <= heterogeneity_threshold):
        
        quality_score = completeness - 5*contamination
        
        filtered_samples.append({
            'Sample': sample_name,
            'Completeness': completeness,
            'Contamination': contamination,
            'Genome_fraction': genome_fraction,
            'N50': n50,
            'Mismatches_per_100kbp': mismatches,
            'Quality_Score': quality_score
        })
        print(f"  {sample_name}: ACCEPTÉ (Score de qualité: {quality_score:.2f})")
    else:
        print(f"  {sample_name}: REJETÉ (C:{completeness:.2f}%, Cont:{contamination:.2f}%, MM:{mismatches:.2f})")

# Créer un DataFrame avec les échantillons filtrés
filtered_df = pd.DataFrame(filtered_samples)

# Si aucun échantillon ne correspond aux critères
if filtered_df.empty:
    print("\nAucun échantillon ne répond aux critères de filtrage.")
    exit(0)

# Trier par score de qualité décroissant
filtered_df = filtered_df.sort_values('Quality_Score', ascending=False)

# Enregistrer les résultats
filtered_df.to_csv('filtered_samples.csv', index=False)

# Afficher les meilleurs échantillons
print("\nTop 5 meilleurs échantillons:")
print(filtered_df.head(5).to_string(index=False))

# Sélectionner le meilleur consensus
best_sample = filtered_df.iloc[0]['Sample']
print(f"\nMeilleur échantillon: {best_sample}")

# Copier le meilleur consensus dans un fichier séparé
shutil.copy(f"consensus/{best_sample}.consensus.fasta", "best_consensus.fasta")
print(f"Le meilleur consensus a été copié dans 'best_consensus.fasta'")
