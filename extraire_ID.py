#!/usr/bin/env python

from Bio import SeqIO

fasta_file = "gisaid_sequences.fasta"
ids_file = "ID_selected.txt"
output_file = "CHIKV_selected.fasta"

# Lire les IDs
with open(ids_file, "r") as f:
    ids_to_keep = set(line.strip() for line in f if line.strip())

print(f"Nombre d'IDs à extraire : {len(ids_to_keep)}")

# Initialiser la liste de résultats
matched_records = []

for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id in ids_to_keep:
        matched_records.append(record)
        print(f"✔️ Séquence trouvée : {record.id}")
    else:
        print(f"❌ Non trouvé : {record.id}")

# Écrire les séquences extraites
if matched_records:
    SeqIO.write(matched_records, output_file, "fasta")
    print(f"✅ {len(matched_records)} séquences extraites dans : {output_file}")
else:
    print("⚠️ Aucune séquence n'a été extraite. Vérifie les noms d'ID.")


