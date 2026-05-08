## en vez de utilizar el jaccard, pruebo con el overlap (https://towardsdatascience.com/similarity-measures-and-graph-adjacency-with-sets-a33d16e527e1/)

library(ggplot2)
library(dplyr)
library(ggrepel)
####
# Function to calculate Jaccard similarity
overlap <- function(x, y) {
  intersect <- sum(x & y) # Count of shared 1s
  min_set <- min(sum(x), sum(y))    # Count of unique 1s
  return(intersect / min_set)
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
  overlap(gene_disease[[i]], gene_disease[[j]])
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

# Compute the number of genes of each disease

n_genes <- colSums(gene_disease==1)
df_n_genes <- as.data.frame(n_genes)
# Compute the number of shared genes for each pair of diseases
n_genes_diff <- matrix(0, ncol = ncol(gene_disease), nrow = ncol(gene_disease))
colnames(n_genes_diff) <- colnames(gene_disease)
rownames(n_genes_diff) <- colnames(gene_disease)
head(n_genes_diff)

for (i in 1:ncol(gene_disease)) {
  for (j in i:ncol(gene_disease)) {
    # Count shared genes with value 1 in both diseases
    diff_count <- sum(gene_disease[, i]) - sum(gene_disease[, j])
    n_genes_diff[i, j] <- diff_count
    n_genes_diff[j, i] <- diff_count # Symmetrical
  }
}

# Convert matrix to data frame for better visualization
df_n_genes_diff <- as.data.frame(n_genes_diff)
head(df_n_genes_diff)
# Print the results
print(df_n_genes_diff)

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
  overlap(gene_disease_glow[[i]], gene_disease_glow[[j]])
}))

#rownames(jaccard_matrix_glow) <- colnames(jaccard_matrix_glow) <- gsub("_GA$", "", colnames(gene_disease_glow))
rownames(jaccard_matrix_glow) <- colnames(jaccard_matrix_glow) <- gene_disease_glow

# Convert to a DataFrame for better readability
jaccard_df_glow <- as.data.frame(jaccard_matrix_glow)
colnames(jaccard_df_glow) <- colnames(gene_disease_glow)
rownames(jaccard_df_glow) <- colnames(gene_disease_glow)

# Print Jaccard Similarity Matrix
print(jaccard_df_glow)

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
  Disease1 = rep(rownames(n_genes_diff), times = ncol(n_genes_diff)),
  Disease2 = rep(colnames(n_genes_diff), each = nrow(n_genes_diff)),
  N_genes_difference = as.vector(n_genes_diff)
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
diagonal_values <- diag(as.matrix(n_genes_diff))

# Get the row names corresponding to the diagonal
row_names <- rownames(n_genes_diff)

merged_results[merged_results$Disease1 == "Intellectual_disability_-_microarray_and_sequencing" ,]$Disease1 <- "Intellectual_disability"
merged_results[merged_results$Disease2 == "Intellectual_disability_-_microarray_and_sequencing" ,]$Disease2 <- "Intellectual_disability"

merged_results[merged_results$Disease1 == "Likely_inborn_error_of_metabolism_-_targeted_testing_not_possible" ,]$Disease1 <- "Likely_inborn_error_of_metabolism"
merged_results[merged_results$Disease2 == "Likely_inborn_error_of_metabolism_-_targeted_testing_not_possible" ,]$Disease2 <- "Likely_inborn_error_of_metabolism"


merged_results <- merged_results %>% dplyr::mutate(my_label = paste(gsub("(_|-)", " " ,Disease1), "--", gsub("_", " " ,Disease2)))

#### vamos a cambiar un poco los tamaños de la diferencia del número de genes porque hay una variedad muy grande (Intellectual disability tiene muchos genes). 
## si la diferencia es > 1000, ponemos 1000 y luego una etiqyeta de >1000. 


merged_results <- merged_results%>% mutate(new_N_genes_diff = case_when((abs(N_genes_difference) >= 1000) ~ 1000,
                                     (abs(N_genes_difference) < 1000 & abs(N_genes_difference) >= 500) ~ 500,
                                     (abs(N_genes_difference) < 500 & abs(N_genes_difference) >= 100) ~ 100,
                                     (abs(N_genes_difference) < 100 & abs(N_genes_difference) >= 50) ~ 50,
                                     (abs(N_genes_difference) < 50 ~ 0)))


subset_data <- rbind(head(merged_results[order(merged_results$Jaccard_difference),],n=5),
                 head(merged_results[order(merged_results$Jaccard_difference, decreasing = TRUE),],n=5))

ggplot(merged_results, aes(x = Jaccard_before_glow, y = Jaccard_after_glow)) +
  geom_point(aes(color = Jaccard_difference, size = abs(N_genes_difference)), alpha = 0.7) +
  geom_text_repel(data = subset_data,  # Specify the subset data
                  aes(x = Jaccard_before_glow, y = Jaccard_after_glow, label = my_label)#,  # Mapping for subset
                  # min.segment.length = 0,  # Ensure no minimum length
                  # segment.size = 0.5,      # Adjust the segment thickness
                  # seed = 42,
                  # max.overlaps = Inf,      # Allow all labels to display
                  # size = 5,
                  # box.padding = unit(1, "lines"),
                  # point.padding = unit(1, "lines")
                  ) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0, name = "Jaccard similarity difference") +
  scale_size_continuous(range = c(1, 8), guide = guide_legend(title = "Gene Set Size Difference"), breaks = c(0,50,250,500,1000,1500)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Scatter Plot of Unique Jaccard Similarities (Before vs After)",
       subtitle = "Colored by Jaccard Difference (After - Before)",
       x = "Before GLOW Similarity",
       y = "After GLOW Similarity") +
  theme_minimal(base_size = 16) + 
  theme(legend.position = "right")




subset_data <- merged_results %>% dplyr::filter(Disease1 == "Confirmed_Fanconi_anaemia_or_Bloom_syndrome" & Disease2 == "Limb_disorders")

ggplot(merged_results, aes(x = Jaccard_before_glow, y = Jaccard_after_glow)) +
  geom_point(aes(color = Jaccard_difference, size = abs(N_genes_difference)), alpha = 0.7) +
  geom_text_repel(data = subset_data,  # Specify the subset data
                  aes(x = Jaccard_before_glow, y = Jaccard_after_glow, label = my_label),
                  min.segment.length = 0
  ) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0, name = "Jaccard similarity difference") +
  scale_size_continuous( guide = guide_legend(title = "Common Genes Difference")) +
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

###################
intersect(glow_predicted_genes[["Mitochondrial_disorder_with_complex_I_deficiency_GA"]], glow_predicted_genes[["Mitochondrial_disorder_with_complex_III_deficiency_GA"]])
intersect(rownames(glowmatrix %>% dplyr::filter(Mitochondrial_disorder_with_complex_I_deficiency_GA == 0)),
          rownames(glowmatrix %>% dplyr::filter(Mitochondrial_disorder_with_complex_III_deficiency_GA == 0)))
###################

###################
intersect(glow_predicted_genes[["Adult_solid_tumours_cancer_susceptibility_GA"]], glow_predicted_genes[["Neuroendocrine_cancer_pertinent_cancer_susceptibility_GA"]])
intersect(rownames(glowmatrix %>% dplyr::filter(Adult_solid_tumours_cancer_susceptibility_GA == 0)),
          rownames(glowmatrix %>% dplyr::filter(Neuroendocrine_cancer_pertinent_cancer_susceptibility_GA == 0)))
###################

length(intersect(glow_predicted_genes[["Confirmed_Fanconi_anaemia_or_Bloom_syndrome_GA"]], glow_predicted_genes[["Limb_disorders_GA"]]))
length(intersect(rownames(glowmatrix %>% dplyr::filter(Confirmed_Fanconi_anaemia_or_Bloom_syndrome_GA == 0)),
          rownames(glowmatrix %>% dplyr::filter(Limb_disorders_GA == 0))))

min(length(glow_predicted_genes[["Confirmed_Fanconi_anaemia_or_Bloom_syndrome_GA"]]), length(glow_predicted_genes[["Limb_disorders_GA"]]))

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
