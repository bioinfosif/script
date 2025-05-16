#!/bin/bash

cd data  # ou le chemin vers ton dossier de FASTQ

for file in *; do
    if [[ "$file" == *-* ]]; then
        new_name="${file//-/0}"
        echo "Renomme : $file → $new_name"
        mv "$file" "$new_name"
    fi
done
