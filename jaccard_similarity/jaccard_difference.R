library(ggplot2)
library(dplyr)
library(ggrepel)
####
# Function to calculate Jaccard similarity
jaccard_similarity <- function(x, y) {
  intersect <- sum(x & y) # Count of shared 1s
  union <- sum(x | y)     # Count of unique 1s
  if (union == 0) return(NA) # Avoid division by zero
  return(intersect / union)
}
####
############## Primero, la matriz de similaridad y el número de genes compartidos en PanelApp:
glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

glowmatrix[glowmatrix>0] <- 1
zero_rows <- apply(glowmatrix, 1, function(row) any(row == 0))
filtered_glow <- glowmatrix[zero_rows, ]
gene_disease <- 1-filtered_glow # Matriz: 1 si el gen está en la lista de PanelApp, 0 si no está

colnames(gene_disease) <- gsub("_GA$", "", colnames(gene_disease)) 

nrow(gene_disease) #4414 # genes únicos en panelapp
ncol(gene_disease) #209 # enfermedades

# Compute Jaccard Similarity Matrix
jaccard_matrix_panelapp <- outer(1:ncol(gene_disease), 1:ncol(gene_disease), Vectorize(function(i, j) {
  jaccard_similarity(gene_disease[[i]], gene_disease[[j]])
}))

#rownames(jaccard_matrix_panelapp) <- colnames(jaccard_matrix_panelapp) <- gsub("_GA$", "", colnames(gene_disease))

rownames(jaccard_matrix_panelapp) <- colnames(jaccard_matrix_panelapp) <- colnames(gene_disease)
# Convert to a DataFrame for better readability
jaccard_df_panelapp <- as.data.frame(jaccard_matrix_panelapp)
colnames(jaccard_df_panelapp) <- colnames(gene_disease)
rownames(jaccard_df_panelapp) <- colnames(gene_disease)

# Print Jaccard Similarity Matrix
print(jaccard_df_panelapp)

## Calculamos el número de genes que comparte cada enfermedad con las demás

# Compute the number of shared genes for each pair of diseases
shared_genes_panelapp <- matrix(0, ncol = ncol(gene_disease), nrow = ncol(gene_disease))
colnames(shared_genes_panelapp) <- colnames(gene_disease)
rownames(shared_genes_panelapp) <- colnames(gene_disease)
head(shared_genes_panelapp)

for (i in 1:ncol(gene_disease)) {
  for (j in i:ncol(gene_disease)) {
    # Count shared genes with value 1 in both diseases
    shared_count <- sum(gene_disease[, i] == 1 & gene_disease[, j] == 1)
    shared_genes_panelapp[i, j] <- shared_count
    shared_genes_panelapp[j, i] <- shared_count # Symmetrical
  }
}


# Convert matrix to data frame for better visualization
shared_genes_df_panelapp <- as.data.frame(shared_genes_panelapp)
head(shared_genes_panelapp)
# Print the results
print(shared_genes_df_panelapp)

####### Ahora, la matriz de similaridad y los genes compartidos tras GLOWgenes

glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

nrow(glowmatrix)
ncol(glowmatrix)

gene_panels <- glowmatrix

gene_panels[gene_panels>0] <- 1
zero_rows <- apply(gene_panels, 1, function(row) any(row == 0))
filtered_glow <- gene_panels[zero_rows, ]
gene_disease <- 1-filtered_glow

nrow(gene_disease) #4414
ncol(gene_disease) #209

gene_panels_glow <- matrix(0, nrow(gene_panels), ncol(gene_panels))
colnames(gene_panels_glow) <- colnames(gene_panels)
rownames(gene_panels_glow) <- rownames(gene_panels)


glowmatrix$genes <- rownames(glowmatrix)

for (disease_i in colnames(gene_panels_glow)){
  n_genes <- length(gene_disease[gene_disease[disease_i] == 1, disease_i])
  genes_glow <- glowmatrix[(glowmatrix[disease_i] <= n_genes & glowmatrix[disease_i] != 0), c("genes", disease_i)] 
  gene_panels_glow[genes_glow$genes, disease_i] <- 1
}

colnames(gene_panels_glow) <- gsub("_GA$", "", colnames(gene_panels_glow)) 

one_rows <- apply(gene_panels_glow, 1, function(row) any(row == 1))
gene_disease_glow <- as.data.frame(gene_panels_glow[one_rows, ])

nrow(gene_disease_glow)
ncol(gene_disease_glow)

# Compute Jaccard Similarity Matrix
jaccard_matrix_glow <- outer(1:ncol(gene_disease_glow), 1:ncol(gene_disease_glow), Vectorize(function(i, j) {
  jaccard_similarity(gene_disease_glow[[i]], gene_disease_glow[[j]])
}))

#rownames(jaccard_matrix_glow) <- colnames(jaccard_matrix_glow) <- gsub("_GA$", "", colnames(gene_disease_glow))
rownames(jaccard_matrix_glow) <- colnames(jaccard_matrix_glow) <- gene_disease_glow

# Convert to a DataFrame for better readability
jaccard_df_glow <- as.data.frame(jaccard_matrix_glow)
colnames(jaccard_df_glow) <- colnames(gene_disease_glow)
rownames(jaccard_df_glow) <- colnames(gene_disease_glow)

# Print Jaccard Similarity Matrix
print(jaccard_df_glow)

# Compute the number of shared genes for each pair of diseases
shared_genes_glow <- matrix(0, ncol = ncol(gene_disease_glow), nrow = ncol(gene_disease_glow))
colnames(shared_genes_glow) <- colnames(gene_disease_glow)
rownames(shared_genes_glow) <- colnames(gene_disease_glow)

for (i in 1:ncol(gene_disease_glow)) {
  for (j in i:ncol(gene_disease_glow)) {
    # Count shared genes with value 1 in both diseases
    shared_count <- sum(gene_disease_glow[, i] == 1 & gene_disease_glow[, j] == 1)
    shared_genes_glow[i, j] <- shared_count
    shared_genes_glow[j, i] <- shared_count # Symmetrical
  }
}

# Convert matrix to data frame for better visualization
shared_genes_glow_df <- as.data.frame(shared_genes_glow)

# Print the results
print(shared_genes_glow_df)

################################## Calculamos la diferencia de similaridad: 


jaccard_difference <- data.frame(
  Disease1 = rep(rownames(jaccard_matrix_panelapp), times = ncol(jaccard_matrix_panelapp)),
  Disease2 = rep(colnames(jaccard_matrix_panelapp), each = nrow(jaccard_matrix_panelapp)),
  Jaccard_before_glow = as.vector(jaccard_matrix_panelapp),
  Jaccard_after_glow = as.vector(jaccard_matrix_glow),
  Jaccard_difference = as.vector(jaccard_matrix_glow - jaccard_matrix_panelapp)
)

unique_jaccard_difference <- jaccard_difference[jaccard_difference$Disease1 < jaccard_difference$Disease2, ]


#View(unique_jaccard_difference)

commongenes_difference <- data.frame(
  Disease1 = rep(rownames(shared_genes_panelapp), times = ncol(shared_genes_panelapp)),
  Disease2 = rep(colnames(shared_genes_panelapp), each = nrow(shared_genes_panelapp)),
  Common_genes_before_glow = as.vector(shared_genes_panelapp),
  Common_genes_after_glow = as.vector(shared_genes_glow),
  Common_genes_difference = as.vector(shared_genes_glow - shared_genes_panelapp)
)

unique_commongenes_difference <- commongenes_difference[commongenes_difference$Disease1 < commongenes_difference$Disease2, ]

merged_results <- merge(unique_jaccard_difference, unique_commongenes_difference, by = c("Disease1", "Disease2"))
#merged_results <- merge(jaccard_difference, commongenes_difference, by = c("Disease1", "Disease2"))

nrow(unique_jaccard_difference)
nrow(unique_commongenes_difference)
nrow(merged_results)

#merged_results <- merged_results %>% dplyr::mutate(my_label = ifelse(Jaccard_difference>0.22 | Jaccard_difference<(-0.4),paste(Disease1, "-", Disease2), ""))

# Add the number of genes of each panel: 

# Extract the diagonal values
diagonal_values <- diag(as.matrix(shared_genes_df_panelapp))

# Get the row names corresponding to the diagonal
row_names <- rownames(shared_genes_df_panelapp)

# Combine the results into a dataframe for clarity
result_n_genes <- data.frame(Disease = row_names, Genes = diagonal_values)

# Print the result
print(result_n_genes) # number of genes of each panel
nrow(merged_results)

colnames(result_n_genes) <- c("Disease", "N_Genes_Disease2")
merged_results <- merge(merged_results, result_n_genes, by.x = "Disease2", by.y = "Disease")

colnames(result_n_genes) <- c("Disease", "N_Genes_Disease1")
merged_results <- merge(merged_results, result_n_genes, by.x = "Disease1", by.y = "Disease")

nrow(merged_results)

merged_results[merged_results$Disease1 == "Intellectual_disability_-_microarray_and_sequencing" ,]$Disease1 <- "Intellectual_disability"
merged_results[merged_results$Disease2 == "Intellectual_disability_-_microarray_and_sequencing" ,]$Disease2 <- "Intellectual_disability"

merged_results[merged_results$Disease1 == "Likely_inborn_error_of_metabolism_-_targeted_testing_not_possible" ,]$Disease1 <- "Likely_inborn_error_of_metabolism"
merged_results[merged_results$Disease2 == "Likely_inborn_error_of_metabolism_-_targeted_testing_not_possible" ,]$Disease2 <- "Likely_inborn_error_of_metabolism"


merged_results <- merged_results %>% dplyr::mutate(my_label = ifelse(Jaccard_difference>0.22 | Jaccard_difference<(-0.4) | abs(Common_genes_difference)>200,
                                                                     paste(gsub("(_|-)", " " ,Disease1), "--", gsub("_", " " ,Disease2)), ""))
#oldK <- GeomLabelRepel$draw_key
# define new key without the text label
#library(grid)
#GeomLabelRepel$draw_key <- function (data, params, size) { draw_key_rect(data) }

# subset_data <- merged_results %>%
#   filter(Jaccard_difference > 0.22 | Jaccard_difference < -0.4 | abs(Common_genes_difference) > 200)

subset_data <- merged_results %>%
  filter(Jaccard_difference > 0.22 | Jaccard_difference < -0.5)


ggplot(merged_results, aes(x = Jaccard_before_glow, y = Jaccard_after_glow)) +
  geom_point(aes(color = Jaccard_difference, size = abs(Common_genes_difference)), alpha = 0.7) +
  geom_text_repel(data = subset_data,  # Specify the subset data
                  aes(x = Jaccard_before_glow, y = Jaccard_after_glow, label = my_label),  # Mapping for subset
                  min.segment.length = 0,  # Ensure no minimum length
                  segment.size = 0.5,      # Adjust the segment thickness
                  seed = 42,  
                  max.overlaps = Inf,      # Allow all labels to display
                  size = 5,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(1, "lines")) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0, name = "Jaccard similarity difference") +
  scale_size_continuous(range = c(1, 10), guide = guide_legend(title = "Common Genes Difference")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Scatter Plot of Unique Jaccard Similarities (Before vs After)",
       subtitle = "Colored by Jaccard Difference (After - Before)",
       x = "Before GLOW Similarity",
       y = "After GLOW Similarity") +
  theme_minimal(base_size = 16) + 
  theme(legend.position = "right")


######## ANÁLISIS DE LOS GENES 

### vamos a sacar los N genes de las predicciones de cada enfermedad

head(glowmatrix)

# Crear una lista vacía para guardar los resultados
glow_predicted_genes <- list()

# Excluyendo la última columna
for (col_name in colnames(glowmatrix)[-ncol(glowmatrix)]) {
  # Obtener la columna actual
  columna <- glowmatrix[[col_name]]
  
  # Contar el número de ceros
  N_ceros <- sum(columna == 0)
  
  # Ordenar la columna
  columna_ordenada <- sort(columna)
  
  # Excluir los valores iguales a cero
  columna_ordenada <- columna_ordenada[columna_ordenada != 0]
  
  # Seleccionar los primeros N_ceros valores
  seleccion <- columna_ordenada[1:N_ceros]
  
  # Obtener los nombres de las filas correspondientes a los valores seleccionados
  nombres_filas <- rownames(glowmatrix)[columna %in% seleccion]
  
  # Guardar el resultado en la lista
  glow_predicted_genes[[col_name]] <- nombres_filas
}
# Imprimir los resultados
print(glow_predicted_genes)


glow_predicted_genes$Bardet_Biedl_syndrome_GA
glow_predicted_genes$Thoracic_dystrophies_GA

### para comparar dos grupos y sacar los genes comunes que ha predicho GLOW:
intersect(glow_predicted_genes[["Bardet_Biedl_syndrome_GA"]], glow_predicted_genes[["Thoracic_dystrophies_GA"]])
union(rownames(glowmatrix %>% dplyr::filter(Bardet_Biedl_syndrome_GA == 0)),
      rownames(glowmatrix %>% dplyr::filter(Thoracic_dystrophies_GA == 0)))

rownames(glowmatrix %>% dplyr::filter(Bardet_Biedl_syndrome_GA == 0))
rownames(glowmatrix %>% dplyr::filter(Thoracic_dystrophies_GA == 0))

#### Ductal_plate_malformation_GA y Polycystic_liver_disease_GA
intersect(glow_predicted_genes[["Ductal_plate_malformation_GA"]], glow_predicted_genes[["Polycystic_liver_disease_GA"]])
glow_predicted_genes[["Ductal_plate_malformation_GA"]]
glow_predicted_genes[["Polycystic_liver_disease_GA"]]

intersect(rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0)),
      rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0)))

rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0))
rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0))

union(rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0)),
          rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0)))

rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0))[! rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0)) %in% intersect(rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0)),
                                                                                             rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0)))]

rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0))[! rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0)) %in% intersect(rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0)),
                                                                                                                                                                     rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0)))]
##########################

rownames(glowmatrix %>% dplyr::filter(Ductal_plate_malformation_GA == 0))
rownames(glowmatrix %>% dplyr::filter(Polycystic_liver_disease_GA == 0))

View(glowmatrix[c("Bardet_Biedl_syndrome_GA","Thoracic_dystrophies_GA")])

intersect(glow_predicted_genes[["Cystic_kidney_disease_GA"]], glow_predicted_genes[["Thoracic_dystrophies_GA"]])
intersect(glow_predicted_genes[["Inherited_predisposition_to_acute_myeloid_leukaemia_(AML)_GA"]], glow_predicted_genes[["Polycystic_liver_disease_GA"]])

glow_predicted_genes[["Inherited_phaeochromocytoma_and_paraganglioma_GA"]]
glow_predicted_genes[["Inherited_phaeochromocytoma_and_paraganglioma_excluding_NF1_GA"]]

intersect(glow_predicted_genes[["Inherited_phaeochromocytoma_and_paraganglioma_GA"]], 
          glow_predicted_genes[["Inherited_phaeochromocytoma_and_paraganglioma_excluding_NF1_GA"]])

union(glow_predicted_genes[["Inherited_phaeochromocytoma_and_paraganglioma_GA"]], 
          glow_predicted_genes[["Inherited_phaeochromocytoma_and_paraganglioma_excluding_NF1_GA"]])

intersect(rownames(glowmatrix %>% dplyr::filter(Inherited_phaeochromocytoma_and_paraganglioma_GA == 0)),
          rownames(glowmatrix %>% dplyr::filter(Inherited_phaeochromocytoma_and_paraganglioma_excluding_NF1_GA == 0)))

union(rownames(glowmatrix %>% dplyr::filter(Inherited_phaeochromocytoma_and_paraganglioma_GA == 0)),
          rownames(glowmatrix %>% dplyr::filter(Inherited_phaeochromocytoma_and_paraganglioma_excluding_NF1_GA == 0)))


#######################################3

panels_classification <- read.delim("~/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/panels/final_reclasified_panels.tsv", header = F, stringsAsFactors = F, quote = "", check.names=F)
colnames(panels_classification) <- c("classification","panel")

unique(panels_classification$classification)
#panels_class <- "Cardiovascular disorders"
panels_class <- "Neurology and neurodevelopmental disorders"


panels_before <- jaccard_df_panelapp[panels_classification[panels_classification$classification == panels_class,]$panel, panels_classification[panels_classification$classification == panels_class,]$panel]

library(pheatmap)

# Convert the dataframe to a matrix
data_matrix <- as.matrix(panels_before)

# Plot the heatmap
pheatmap(data_matrix,
         cluster_rows = TRUE,   # Clustering rows
         cluster_cols = TRUE,   # Clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(50), # Color scale
         main = "Heatmap of Conditions", # Title
         fontsize_row = 8,      # Font size for row labels
         fontsize_col = 8)      # Font size for column labels



panels_after <- jaccard_df_glow[panels_classification[panels_classification$classification == panels_class,]$panel, panels_classification[panels_classification$classification == panels_class,]$panel]

data_matrix_after <- as.matrix(panels_after)

# Plot the heatmap
pheatmap(data_matrix_after,
         cluster_rows = TRUE,   # Clustering rows
         cluster_cols = TRUE,   # Clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(50), # Color scale
         main = "Heatmap of Conditions", # Title
         fontsize_row = 8,      # Font size for row labels
         fontsize_col = 8)      # Font size for column labels

## no clustering


pheatmap(data_matrix,
         cluster_rows = FALSE,   # Clustering rows
         cluster_cols = FALSE,   # Clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(50), # Color scale
         main = "Heatmap of Conditions", # Title
         fontsize_row = 8,      # Font size for row labels
         fontsize_col = 8)      # Font size for column labels

pheatmap(data_matrix_after,
         cluster_rows = FALSE,   # Clustering rows
         cluster_cols = FALSE,   # Clustering columns
         color = colorRampPalette(c("blue", "white", "red"))(50), # Color scale
         main = "Heatmap of Conditions", # Title
         fontsize_row = 8,      # Font size for row labels
         fontsize_col = 8)      # Font size for column labels

