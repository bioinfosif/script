# Chargement des packages nécessaires
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# Lecture du fichier TSV
abundance <- read.table("bracken_abundance_table.tsv", header=TRUE, sep="\t", row.names=1)

# Diagnostic des données
print("Dimensions des données:")
print(dim(abundance))
print("\nStructure des données:")
str(abundance)
print("\nPremières lignes:")
head(abundance)

# Vérification des valeurs non-numériques
print("\nVérification des colonnes non-numériques:")
non_numeric_cols <- sapply(abundance, function(x) !is.numeric(x))
print(names(abundance)[non_numeric_cols])

# Conversion en matrice numérique avec gestion des erreurs
convert_to_numeric <- function(x) {
  # Remplacer les valeurs vides, NA, ou non-numériques par 0
  x[is.na(x) | x == "" | x == " "] <- 0
  as.numeric(as.character(x))
}

# Appliquer la conversion à toutes les colonnes
abundance_numeric <- as.data.frame(lapply(abundance, convert_to_numeric))
rownames(abundance_numeric) <- rownames(abundance)

# Remplacer les NA par 0 après conversion
abundance_numeric[is.na(abundance_numeric)] <- 0

# Vérification après conversion
print("\nAprès conversion numérique:")
print(paste("Nombre de NA:", sum(is.na(abundance_numeric))))
print(paste("Nombre de valeurs infinies:", sum(is.infinite(as.matrix(abundance_numeric)))))

# Filtrage des espèces avec au moins une valeur > 0.01 dans au moins un échantillon
filtered_abundance <- abundance_numeric[rowSums(abundance_numeric > 0.01, na.rm = TRUE) > 0, ]

print(paste("\nNombre d'espèces après filtrage:", nrow(filtered_abundance)))
print(paste("Nombre d'échantillons:", ncol(filtered_abundance)))

# Si aucune espèce ne passe le filtre 0.01, essayer avec un seuil plus bas
if(nrow(filtered_abundance) == 0) {
  print("Aucune espèce avec abondance > 0.01. Essai avec seuil 0.001...")
  filtered_abundance <- abundance_numeric[rowSums(abundance_numeric > 0.001, na.rm = TRUE) > 0, ]
  print(paste("Nombre d'espèces avec seuil 0.001:", nrow(filtered_abundance)))
}

# Si encore aucune espèce, prendre les top espèces par somme
if(nrow(filtered_abundance) == 0) {
  print("Prise des 20 espèces les plus abondantes...")
  row_sums <- rowSums(abundance_numeric, na.rm = TRUE)
  top_species <- names(sort(row_sums, decreasing = TRUE))[1:min(20, length(row_sums))]
  filtered_abundance <- abundance_numeric[top_species, ]
}

# Conversion en matrice pour pheatmap
mat <- as.matrix(filtered_abundance)

# Vérifications finales avant heatmap
print("\nVérifications finales:")
print(paste("Dimensions de la matrice:", paste(dim(mat), collapse = " x ")))
print(paste("Valeurs min/max:", min(mat, na.rm = TRUE), "/", max(mat, na.rm = TRUE)))

# Supprimer les lignes avec toutes les valeurs à 0
mat <- mat[rowSums(mat, na.rm = TRUE) > 0, ]

# Supprimer les colonnes avec toutes les valeurs à 0
mat <- mat[, colSums(mat, na.rm = TRUE) > 0]

print(paste("Après suppression des lignes/colonnes nulles:", paste(dim(mat), collapse = " x ")))

# Fonction pour créer le heatmap avec gestion d'erreurs
create_heatmap <- function(matrix_data, scale_type = "row") {
  tryCatch({
    # Vérifier si on peut faire le clustering
    if(nrow(matrix_data) > 1 && ncol(matrix_data) > 1) {
      pheatmap(
        matrix_data,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        scale = scale_type,
        color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
        fontsize_row = 8,
        fontsize_col = 10,
        main = "Heatmap des espèces (Bracken)",
        na_col = "white",
        fontface = "italic"
      )
    } else {
      # Si trop peu de données pour clustering
      pheatmap(
        matrix_data,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        scale = scale_type,
        color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
        fontsize_row = 8,
        fontsize_col = 10,
        main = "Heatmap des espèces (Bracken) - Sans clustering",
        na_col = "white",
        fontface = "italic"
      )
    }
  }, error = function(e) {
    print(paste("Erreur avec scale =", scale_type, ":", e$message))
    
    # Essayer sans scaling
    if(scale_type != "none") {
      print("Tentative sans scaling...")
      create_heatmap(matrix_data, "none")
    } else {
      print("Impossible de créer le heatmap. Affichage des données de base:")
      print(head(matrix_data))
    }
  })
}

# Créer le heatmap
print("\nCréation du heatmap...")
create_heatmap(mat, "row")

# Alternative: heatmap simple si pheatmap ne fonctionne pas
print("\nAlternative avec heatmap() de base:")
tryCatch({
  heatmap(mat, 
          col = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
          main = "Heatmap alternatif des espèces (Bracken)")
}, error = function(e) {
  print(paste("Erreur avec heatmap() aussi:", e$message))
})

# Affichage des statistiques finales
print("\nStatistiques finales:")
print(paste("Nombre total d'espèces dans la matrice finale:", nrow(mat)))
print(paste("Nombre total d'échantillons:", ncol(mat)))
print("Espèces incluses:")
print(head(rownames(mat), 10))
