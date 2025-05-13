#!/usr/bin/env python


# Dictionnaire de correspondance
id_mapping = {
    "2": "IRESSEF_CHIK0048/2023",
    "3": "IRESSEF_CHIK0056/2023",
    "9": "IRESSEF_CHIK0092/2023",
    "11": "IRESSEF_CHIK0120/2023",
    "12": "IRESSEF_CHIK0134/2023",
    "13": "IRESSEF_CHIK0136/2023",
    "17": "IRESSEF_CHIK0167/2023",
    "18": "IRESSEF_CHIK0171/2023",
    "19": "IRESSEF_CHIK0179/2023",
    "20": "IRESSEF_CHIK0191/2023",
    "22": "IRESSEF_CHIK0210/2023",
    "23": "IRESSEF_CHIK0219/2023",
    "25": "IRESSEF_CHIK0238/2023",
    "26": "IRESSEF_CHIK0241/2023",
    "27": "IRESSEF_CHIK0250/2023",
    "28": "IRESSEF_CHIK0252/2023",
    "29": "IRESSEF_CHIK0286/2023",
    "31": "IRESSEF_CHIK0307/2023",
    "33": "IRESSEF_CHIK0333/2023",
    "34": "IRESSEF_CHIK0339/2023",
    "35": "IRESSEF_CHIK0345/2023",
    "37": "IRESSEF_CHIK0358/2023",
    "38": "IRESSEF_CHIK0372/2023",
    "39": "IRESSEF_CHIK0386/2023",
    "40": "IRESSEF_CHIK0398/2023",
    "42": "IRESSEF_CHIK0424/2023",
    "43": "IRESSEF_CHIK0425/2023",
    "45": "IRESSEF_CHIK0443/2023",
    "47": "IRESSEF_CHIK0470/2023",
    "50": "IRESSEF_CHIK0499/2023",
    "51": "IRESSEF_CHIK0503/2023"
}

# Fichiers d'entrée et de sortie
input_fasta = "CHIKV_consensus.fasta"      # Remplace par ton fichier original
output_fasta = "renamed.fasta"   # Fichier renommé

# Lecture et écriture avec remplacement
with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            original_id = line[1:].strip()
            new_id = id_mapping.get(original_id, original_id)
            outfile.write(f">{new_id}\n")
        else:
            outfile.write(line)

