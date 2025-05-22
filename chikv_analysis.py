#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
from collections import defaultdict

# Dossier contenant vos fichiers .tabular
dossier = "./Report"  # Ajustez selon votre structure

# Fonction pour extraire le nom de l'échantillon
def get_sample_name(filename):
    # Format attendu: NNN_SNN_L001_RN_001.fastq.gz.tabular
    match = re.match(r'(.+_S\d+)_L001_R[12]_001\.fastq\.gz\.tabular', filename)
    if match:
        return match.group(1)  # Retourne la partie avant "_L001_R"
    return None

# Fonction pour trouver le pourcentage de CHIKV dans un fichier
def find_chikv_percentage(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            # Recherche de différentes appellations possibles du virus Chikungunya
            if "Chikungunya virus" in line or "CHIKV" in line or "alphavirus chikungunya" in line:
                elements = line.strip().split('\t')
                if len(elements) >= 2:
                    try:
                        # Le pourcentage est dans la première colonne
                        percentage = float(elements[0])
                        return percentage
                    except ValueError:
                        pass
    return 0.0  # Retourne 0 si aucune entrée CHIKV n'est trouvée

# Regrouper les fichiers par échantillon
samples = defaultdict(list)
for filename in os.listdir(dossier):
    if filename.endswith(".tabular"):
        sample_name = get_sample_name(filename)
        if sample_name:
            samples[sample_name].append(os.path.join(dossier, filename))

# Analyser chaque échantillon
results = []
for sample_name, files in samples.items():
    chikv_percentages = []
    
    for file in files:
        percentage = find_chikv_percentage(file)
        if percentage > 0:
            chikv_percentages.append(percentage)
            # Afficher les détails pour déboguer
            print(f"Fichier: {os.path.basename(file)}, CHIKV: {percentage:.2f}%")
    
    # Calculer la moyenne des proportions entre R1 et R2
    if chikv_percentages:
        avg_percentage = sum(chikv_percentages) / len(chikv_percentages)
        results.append((sample_name, avg_percentage, chikv_percentages))
    else:
        # Si aucun CHIKV détecté, ajouter avec 0%
        results.append((sample_name, 0.0, []))

# Trier les résultats par proportion décroissante
results.sort(key=lambda x: x[1], reverse=True)

# Afficher tous les résultats
print("\nRésultats pour tous les échantillons:")
print("-" * 60)
print(f"{'Échantillon':<15} | {'% CHIKV (moyen)':<15} | {'% CHIKV (R1, R2)':<20}")
print("-" * 60)
for sample, avg_percentage, percentages in results:
    percentages_str = ", ".join([f"{p:.2f}%" for p in percentages]) if percentages else "0%"
    print(f"{sample:<15} | {avg_percentage:>12.2f}% | {percentages_str:<20}")

# Afficher uniquement les échantillons avec plus de 15% de CHIKV
above_threshold = [r for r in results if r[1] > 15]
if above_threshold:
    print("\nÉchantillons avec plus de 15% de CHIKV:")
    print("-" * 60)
    print(f"{'Échantillon':<15} | {'% CHIKV (moyen)':<15} | {'% CHIKV (R1, R2)':<20}")
    print("-" * 60)
    for sample, avg_percentage, percentages in above_threshold:
        percentages_str = ", ".join([f"{p:.2f}%" for p in percentages])
        print(f"{sample:<15} | {avg_percentage:>12.2f}% | {percentages_str:<20}")
else:
    print("\nAucun échantillon ne contient plus de 15% de CHIKV.")
