from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Charger l'alignement
alignment = AlignIO.read("Final_All_CHIKV_aligned_trim.fasta", "fasta")

# Obtenir la longueur de l'alignement
align_len = alignment.get_alignment_length()

# Trouver les positions variables
variable_positions = []
for i in range(align_len):
    column = alignment[:, i]
    unique = set(column.replace("-", ""))  # Ignorer les gaps
    if len(unique) > 1:
        variable_positions.append(i)

# Construire un nouvel alignement avec uniquement les positions variables
new_records = []
for record in alignment:
    new_seq = ''.join([record.seq[i] for i in variable_positions])
    record.seq = new_seq
    new_records.append(record)

# Sauvegarder l’alignement nettoyé
AlignIO.write(MultipleSeqAlignment(new_records), "alignment_variable_only.fasta", "fasta")

print(f"{align_len - len(variable_positions)} colonnes conservées supprimées.")
