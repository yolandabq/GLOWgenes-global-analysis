color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
colnames(color_panels)<-c("type_panel","panel")
color_panels$type_panel<-ifelse(color_panels$type_panel %in% c("Growth disorders",
                                                               "Hearing and ear disorders","Rheumatological disorders",
                                                               "Viral research"), "Other_panels",color_panels$type_panel)

glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

glowmatrix[glowmatrix>0] <- 1
zero_rows <- apply(glowmatrix, 1, function(row) any(row == 0))
filtered_glow <- glowmatrix[zero_rows, ]
gene_disease <- 1-filtered_glow

nrow(gene_disease) #4414
ncol(gene_disease) #209



# Function to calculate Jaccard similarity
jaccard_similarity <- function(x, y) {
  intersect <- sum(x & y) # Count of shared 1s
  union <- sum(x | y)     # Count of unique 1s
  if (union == 0) return(NA) # Avoid division by zero
  return(intersect / union)
}

# Compute Jaccard Similarity Matrix
jaccard_matrix <- outer(1:ncol(gene_disease), 1:ncol(gene_disease), Vectorize(function(i, j) {
  jaccard_similarity(gene_disease[[i]], gene_disease[[j]])
}))

# Convert to a DataFrame for better readability
jaccard_df <- as.data.frame(jaccard_matrix)
colnames(jaccard_df) <- colnames(gene_disease)
rownames(jaccard_df) <- colnames(gene_disease)

# Print Jaccard Similarity Matrix
print(jaccard_df)

library(igraph)

rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- gsub("_GA$", "", colnames(gene_disease))
jaccard_matrix_diseases <- jaccard_matrix

# Compute the number of shared genes for each pair of diseases
shared_genes <- matrix(0, ncol = ncol(gene_disease), nrow = ncol(gene_disease))
colnames(shared_genes) <- colnames(gene_disease)
rownames(shared_genes) <- colnames(gene_disease)

for (i in 1:ncol(gene_disease)) {
  for (j in i:ncol(gene_disease)) {
    # Count shared genes with value 1 in both diseases
    shared_count <- sum(gene_disease[, i] == 1 & gene_disease[, j] == 1)
    shared_genes[i, j] <- shared_count
    shared_genes[j, i] <- shared_count # Symmetrical
  }
}

# Convert matrix to data frame for better visualization
shared_genes_diseases_df <- as.data.frame(shared_genes)

# Print the results
print(shared_genes_diseases_df)

# Convert the Jaccard matrix to a data frame of edges
edges <- which(jaccard_matrix > 0, arr.ind = TRUE) # Find all non-zero entries
edges_df <- data.frame(
  from = rownames(jaccard_matrix)[edges[, 1]],
  to = colnames(jaccard_matrix)[edges[, 2]],
  weight = jaccard_matrix[edges] # Jaccard Index as edge weight
)

# Optional: Remove self-loops (CARDIOT Y CARDIOT QUE TIENEN UN INDICE DE 1) and duplicate edges (since it's symmetric) 
#-> cardiot a benign cat es igual a benign cat a cardiot
edges_df <- edges_df[edges_df$from != edges_df$to, ]
edges_df <- edges_df[!duplicated(t(apply(edges_df, 1, sort))), ]

# add interaction type = ddi (disease - disease to interaction)
edges_df$interaction<-"ddi"
#poner interaction en el medio
edges_df<-edges_df[,c("from","interaction","to","weight")]

#write.table(edges_df,file="~/tblab/yolanda/GLOWgenes/panelAPP/analysis/jaccard_network/jaccard_diseases_edges.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

###############3

### comprobar todo esto bien


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

nrow(gene_panels_glow)
ncol(gene_panels_glow)

one_rows <- apply(gene_panels_glow, 1, function(row) any(row == 1))
gene_disease <- as.data.frame(gene_panels_glow[one_rows, ])

nrow(gene_disease)
ncol(gene_disease)

# Compute Jaccard Similarity Matrix
jaccard_matrix <- outer(1:ncol(gene_disease), 1:ncol(gene_disease), Vectorize(function(i, j) {
  jaccard_similarity(gene_disease[[i]], gene_disease[[j]])
}))

# Convert to a DataFrame for better readability
jaccard_df <- as.data.frame(jaccard_matrix)
colnames(jaccard_df) <- colnames(gene_disease)
rownames(jaccard_df) <- colnames(gene_disease)

# Print Jaccard Similarity Matrix
print(jaccard_df)

library(igraph)

rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- gsub("_GA$", "", colnames(gene_disease))
jaccard_matrix_glow <- jaccard_matrix

# Compute the number of shared genes for each pair of diseases
shared_genes <- matrix(0, ncol = ncol(gene_disease), nrow = ncol(gene_disease))
colnames(shared_genes) <- colnames(gene_disease)
rownames(shared_genes) <- colnames(gene_disease)

for (i in 1:ncol(gene_disease)) {
  for (j in i:ncol(gene_disease)) {
    # Count shared genes with value 1 in both diseases
    shared_count <- sum(gene_disease[, i] == 1 & gene_disease[, j] == 1)
    shared_genes[i, j] <- shared_count
    shared_genes[j, i] <- shared_count # Symmetrical
  }
}

# Convert matrix to data frame for better visualization
shared_genes_glow_df <- as.data.frame(shared_genes)

# Print the results
print(shared_genes_glow_df)

# Convert the Jaccard matrix to a data frame of edges
edges <- which(jaccard_matrix > 0, arr.ind = TRUE) # Find all non-zero entries
edges_df <- data.frame(
  from = rownames(jaccard_matrix)[edges[, 1]],
  to = colnames(jaccard_matrix)[edges[, 2]],
  weight = jaccard_matrix[edges] # Jaccard Index as edge weight
)

# Optional: Remove self-loops (CARDIOT Y CARDIOT QUE TIENEN UN INDICE DE 1) and duplicate edges (since it's symmetric) 
#-> cardiot a benign cat es igual a benign cat a cardiot
edges_df <- edges_df[edges_df$from != edges_df$to, ]
edges_df <- edges_df[!duplicated(t(apply(edges_df, 1, sort))), ]

# add interaction type = ddi (disease - disease to interaction)
edges_df$interaction<-"ddi"
#poner interaction en el medio
edges_df<-edges_df[,c("from","interaction","to","weight")]


#write.table(edges_df,file="~/tblab/yolanda/GLOWgenes/panelAPP/analysis/jaccard_network/jaccard_glow_edges.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

##################################

difference <- jaccard_matrix_glow - jaccard_matrix_diseases

indices <- which(difference > 0, arr.ind = TRUE)

# Extract results
results <- data.frame(
  Disease1 = rownames(difference)[indices[, 1]],
  Disease2 = colnames(difference)[indices[, 2]],
  Before = apply(indices, 1, function(idx) jaccard_matrix_diseases[idx[1], idx[2]]),
  After = apply(indices, 1, function(idx) jaccard_matrix_glow[idx[1], idx[2]]),
  Increase = apply(indices, 1, function(idx) difference[idx[1], idx[2]])
)

View(results)


ggplot(results, aes(x = Before, y = After)) +
  geom_point(aes(color = Increase > 0), size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Comparison of Jaccard Similarities (Before vs After)",
       x = "Before Similarity",
       y = "After Similarity") +
  theme_minimal()




all_differences <- data.frame(
  Disease1 = rep(rownames(jaccard_matrix_diseases), times = ncol(jaccard_matrix_diseases)),
  Disease2 = rep(colnames(jaccard_matrix_diseases), each = nrow(jaccard_matrix_diseases)),
  Before = as.vector(jaccard_matrix_diseases),
  After = as.vector(jaccard_matrix_glow),
  Difference = as.vector(jaccard_matrix_glow - jaccard_matrix_diseases)
)


# Scatter plot of all differences
ggplot(all_differences, aes(x = Before, y = After)) +
  geom_point(aes(color = Difference, size = abs(Difference)), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 5), guide = guide_legend(title = "Abs Difference")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Scatter Plot of Jaccard Similarities (Before vs After)",
       subtitle = "Colored by Difference (After - Before)",
       x = "Before Similarity",
       y = "After Similarity") +
  theme_minimal() +
  theme(legend.position = "right")


# Remove duplicates by keeping only unique pairs (Disease1 < Disease2)
unique_differences <- all_differences[all_differences$Disease1 < all_differences$Disease2, ]

View(unique_differences)
# Scatter plot of all differences (unique pairs)
ggplot(unique_differences, aes(x = Before, y = After)) +
  geom_point(aes(color = Difference, size = abs(Difference)), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 5), guide = guide_legend(title = "Abs Difference")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Scatter Plot of Unique Jaccard Similarities (Before vs After)",
       subtitle = "Colored by Difference (After - Before)",
       x = "Before Similarity",
       y = "After Similarity") +
  theme_minimal() +
  theme(legend.position = "right")




shared_genes_diseases_df["Cystic_kidney_disease_GA", "Thoracic_dystrophies_GA"]
shared_genes_glow_df["Cystic_kidney_disease_GA", "Thoracic_dystrophies_GA"]

shared_genes_diseases_df["Bardet_Biedl_syndrome_GA", "Thoracic_dystrophies_GA"]
shared_genes_glow_df["Bardet_Biedl_syndrome_GA", "Thoracic_dystrophies_GA"]


shared_genes_diseases_df["Bardet_Biedl_syndrome_GA", "Cystic_kidney_disease_GA"]
shared_genes_glow_df["Bardet_Biedl_syndrome_GA", "Cystic_kidney_disease_GA"]


difference <- shared_genes_glow_df - shared_genes_diseases_df
rownames(shared_genes_glow_df) <- colnames(shared_genes_glow_df) <- gsub("_GA$", "", colnames(shared_genes_glow_df))
rownames(shared_genes_diseases_df) <- colnames(shared_genes_diseases_df) <- gsub("_GA$", "", colnames(shared_genes_diseases_df))

shared_genes_glow <- as.matrix(shared_genes_glow_df)
shared_genes_diseases <- as.matrix(shared_genes_diseases_df)

results_shared_genes <- data.frame(
  Disease1 = rep(rownames(shared_genes_diseases), times = ncol(shared_genes_diseases)),
  Disease2 = rep(colnames(shared_genes_diseases), each = nrow(shared_genes_diseases)),
  Before = as.vector(shared_genes_diseases),
  After = as.vector(shared_genes_glow),
  Difference = as.vector(shared_genes_glow - shared_genes_diseases)
)


unique_shared_genes_differences <- results_shared_genes[results_shared_genes$Disease1 < results_shared_genes$Disease2, ]

merged_results <- merge(all_differences, results_shared_genes, by = c("Disease1", "Disease2"))
unique_shared_genes_differences <- merged_results[merged_results$Disease1 < merged_results$Disease2, ]

ggplot(unique_shared_genes_differences, aes(x = Before.x, y = After.x)) +
  geom_point(aes(color = Difference.x, size = abs(Difference.y)), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(1, 5), guide = guide_legend(title = "Abs Difference")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Scatter Plot of Unique Jaccard Similarities (Before vs After)",
       subtitle = "Colored by Difference (After - Before)",
       x = "Before Similarity",
       y = "After Similarity") +
  theme_minimal() +
  theme(legend.position = "right")




################## cargamos el clustering que nos ha hecho Cytoscape

cyto_cluster_diseases <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/jaccard_network/jaccard_diseases_clustering.csv", header = T, check.names=F, sep = ",")
cyto_cluster_glow <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/jaccard_network/jaccard_glow_clustering.csv", header = T, check.names=F, sep = ",")

cyto_cluster_glow %>% dplyr::filter(`__mclCluster`== 1)
cyto_cluster_diseases %>% dplyr::filter(`__mclCluster`== 1)

## hacemos el merge? y vemos si comparten los clusters de la misma forma 





