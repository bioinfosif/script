#!/bin/bash

CONSENSUS_DIR="CHIKV_out/consensus"

for file in "$CONSENSUS_DIR"/*_consensus.fasta; do
    filename=$(basename "$file")
    new_header=$(echo "$filename" | cut -d'_' -f1)
    sed -i "1s/.*/>${new_header}/" "$file"
    echo "✅ Entête modifié pour $filename --> >$new_header"
done

echo "🎉 Tous les entêtes ont été mis à jour."
