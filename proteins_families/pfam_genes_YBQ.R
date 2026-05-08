
library(utils)
library(R.utils)
library(vcfR)
library(stats)
library(dplyr)
library(readxl)

library(dplyr)
library(factoextra)
library(FactoMineR)
library(scales)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(scales)
library(randomcoloR)
library(ggplot2)
library(ggrepel)
library(gridExtra)

###READ MATRIX DE GENES POR PANELES
glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#extraer lista de genes: hgnc
all_genes<-rownames(glowmatrix)
length(unique(all_genes)) #verificar que son unicos: 24757

########################### 4) FINALMENTE INVESTIGUE PANTHER Y ME QUEDE CON ELLA ##################

################## PANTHER: https://pantherdb.org/
## 1) BUSCAR LOS DATOS: fue bastante lio porque no los encontraba, al final desde la web: https://pantherdb.org/
## me fui a browse>species, en species seleccione homo sapiens > y le di al cuadrado de 20580 genes que hay a la derecha en azul
## Despues sale una tabla de 5 columnas, borre la de "species", me quede con 4 columnas y le di a display: 25000 para que me salieran
## todos los genes, y ya con eso le di a send list > file y ese es el archivo que me descargue.
# HA CAMBIADO UN PELIN ALGO DE INFO PERO LAS PROTEIN CLASSES SON IGUALES

#https://pantherdb.org/list/list.do?filterLevel=1&sortField=PANTHER_PROTEIN_CLASS&listType=1&sortOrder=1&trackingId=5F16A574EC3D283B806C127851D00516&species=All&save=yes&basketItems=all
setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")

############## preparar archivo panther_gene_list ############
original_panther_gene_list <- read.delim("pantherGeneList.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) #All 20589 human genes
panther_gene_list <- read.delim("pantherGeneList.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) #All 20589 human genes
colnames(panther_gene_list)<-c("GeneID","Various","PANTHERFamilySubfamilyIDName","PANTHERProteinClass")
panther_gene_list <- tidyr::separate(panther_gene_list, "Various", into = c("GeneName", "GeneSymbol", "PersistentID","Orthologs"), sep = ";")

# Extracting parts before and after "("
panther_gene_list$PANTHERSubfamilyName <- sub(" \\(.*", "", panther_gene_list$PANTHERFamilySubfamilyIDName) #separe FamilySUbfamilyIDName
panther_gene_list$PANTHERFamilySubfamilyID <- sub(".*\\(", "(", panther_gene_list$PANTHERFamilySubfamilyIDName) #separe FamilySUbfamilyIDName
panther_gene_list$PANTHERFamilySubfamilyID<- gsub("[()]", "", panther_gene_list$PANTHERFamilySubfamilyID) # Remove "(" and ")"

panther_gene_list <- panther_gene_list[, -c(6)] #remove original column: PANTHERFamilySubfamilyIDName
panther_gene_list <- panther_gene_list[, c(1:5, 7:8, 6)] #reorder

panther_gene_list$PANTHERProteinClass<- gsub("\\(", " (", panther_gene_list$PANTHERProteinClass) #add space before "()
panther_gene_list$PANTHERProteinClassName_original <- sub(" \\(.*", "", panther_gene_list$PANTHERProteinClass) #separe FamilySUbfamilyIDName
panther_gene_list$PANTHERProteinClassID_original <- sub(".*\\(", "(", panther_gene_list$PANTHERProteinClass)
panther_gene_list$PANTHERProteinClassID_original<- gsub("[()]", "", panther_gene_list$PANTHERProteinClassID_original) # Remove "(" and ")"
panther_gene_list <- panther_gene_list[, -c(8)] #remove original column: PANTHERProteinClass


############## preparar archivo relationship_protein_class ###############
###  PC00000 es la protein class basica que contiene 34 subclases (a la derecha la )
###http://data.pantherdb.org/PANTHER19.0/ontology/Protein_class_relationship
relationship_protein_class <- read.delim("protein_class_relationship.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) #family and genes
relationship_protein_class<-relationship_protein_class[-c(1,2),]
colnames(relationship_protein_class)<-c("RelationProteinClassID_level2","RelationProteinClassName_level2","RelationProteinClassID_level1","RelationProteinClassName_level1","IDofLevel2inLevel1")
unique_relationship_protein_class<-length(unique(relationship_protein_class$V4)) #6971 familias de proteinas


#### comprobar que todas las: PANTHERProteinClassID_original (archivo: panther_gene_list) están en las "branches" de la relationship (lado de la izquierda del relationship)
unique_pc_genelist<-unique(panther_gene_list$PANTHERProteinClassID_original)
unique_pc_relationship<-unique(relationship_protein_class$RelationProteinClassID_level2)
unique_pc_genelist %in% unique_pc_relationship

#si machean los branches "nivel 2 de class protein" entonces mutate el nivel 1 (la derecha)
merged_df <- merge(panther_gene_list, relationship_protein_class, by.x = "PANTHERProteinClassID_original", by.y = "RelationProteinClassID_level2", all.x = TRUE)
merged_df<-merged_df[,c(2:9,1,10:13)]
merged_df<-merged_df[,-c(13)]
names(merged_df)[names(merged_df) == 'RelationProteinClassName_level2'] <- 'level3_name'
names(merged_df)[names(merged_df) == 'RelationProteinClassName_level1'] <- 'level2_name'
names(merged_df)[names(merged_df) == 'PANTHERProteinClassID_original'] <- 'level3_ID'
names(merged_df)[names(merged_df) == 'RelationProteinClassID_level1'] <- 'level2_ID'

#ordenar dataframe
merged_df <- merged_df[order(merged_df$GeneSymbol,na.last = TRUE),]

#volver a hacer el match porque algunas original estaban más brancheadas 
merged_df_new <- merge(merged_df, relationship_protein_class, by.x = "level2_ID", by.y = "RelationProteinClassID_level2",all.x = TRUE)
merged_df_new <- merged_df_new[order(merged_df_new$GeneSymbol,na.last = TRUE),]
merged_df_new<-merged_df_new[,c(2:11,1,12:16)]
merged_df_new<-merged_df_new[,-c(16)]
merged_df_new<-merged_df_new[,-c(13)] #quitar columna repe al haber hecho el match
names(merged_df_new)[names(merged_df_new) == 'RelationProteinClassID_level1'] <- 'level1_ID'
names(merged_df_new)[names(merged_df_new) == 'RelationProteinClassName_level1'] <- 'level1_name'

#subir el protein class del nivel 2 al nivel 1 (son NAs) -> convertir los niveles 1
merged_df_new$level1_name[is.na(merged_df_new$level1_name)] <- merged_df_new$level2_name[is.na(merged_df_new$level1_name)]
merged_df_new$level1_ID[is.na(merged_df_new$level1_ID)] <- merged_df_new$level2_ID[is.na(merged_df_new$level1_ID)]
### convertir los niveles 2 que son protein class (subirles el nivel 3)
merged_df_new$level2_name[merged_df_new$level2_name == "protein class"] <- merged_df_new$level3_name[merged_df_new$level2_name== "protein class"]
merged_df_new$level2_ID[merged_df_new$level2_ID == "PC00000"] <- merged_df_new$level3_ID[merged_df_new$level2_ID== "PC00000"]

### convertir los niveles 1 que son protein class (subirles el nivel 2)
merged_df_new$level1_name[merged_df_new$level1_name == "protein class"] <- merged_df_new$level2_name[merged_df_new$level1_name== "protein class"]
merged_df_new$level1_ID[merged_df_new$level1_ID == "PC00000"] <- merged_df_new$level2_ID[merged_df_new$level1_ID== "PC00000"]

###si el level 1 machea con el level 2 entonces al level 2 subele el level 3 y el level 3 ponerle -----
non_empty <- merged_df_new$PANTHERProteinClassName_original != ""
same_values <- merged_df_new$level2_ID== merged_df_new$level1_ID
merged_df_new$level2_ID[non_empty & same_values] <- merged_df_new$level3_ID[same_values & non_empty]
merged_df_new$level3_ID[non_empty & same_values] <- "-------"

merged_df_new$level2_name[non_empty & same_values] <- merged_df_new$level3_name[same_values & non_empty]
merged_df_new$level3_name[non_empty & same_values] <- "-------"

############### si el level 2 machea con el level 1 entonces el level 2 ponerle rayitas ####################

same_values <- merged_df_new$level2_ID== merged_df_new$level1_ID
merged_df_new$level2_name[non_empty & same_values] <- "--------"
merged_df_new$level2_ID[non_empty & same_values] <- "--------"

####### a los genes sin info de familia ponerles NO FAMILY
length(unique(merged_df_new$level1_name))
merged_df_new$level1_name <- ifelse(merged_df_new$PANTHERProteinClassName_original== "","NO FAMILY",merged_df_new$level1_name)
table(merged_df_new$level1_name) #6691 genes sin info de familia
merged_df_new
#write.table(merged_df_new, file = "panther_gene_protein_class_gur.txt", sep = "\t", row.names = FALSE)
############################# FIN ###################################


########################################### PCA GLOW MATRIX #####################################
setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")
library(dplyr)
library(factoextra)
library(FactoMineR)
library(scales)
library(randomcoloR)
library(ggplot2)
library(ggrepel)
library(gridExtra)

#1) Read the tab-separated file back into R as a dataframe
gur_parsed_df <- read.table("panther_gene_protein_class_gur.txt", header = TRUE, sep = "\t")
gur_parsed_df <- gur_parsed_df %>%
  dplyr::filter(!GeneSymbol=="unassigned") #20364/20589
##1a) PROTEIN LEVEL 1
genes_class<-gur_parsed_df[,c(3,14)] #protein level 1
genes_class <- unique(genes_class)
genes_class<-genes_class[-c(6835,10211),] #GPX6,MKKS remove el no family que es una proteina desconocida con el mismo nomnbre de lo original
rownames(genes_class)<-genes_class[,1]
genes_class <- genes_class[, c("level1_name"), drop = FALSE]

## 1b) PROTEIN LEVEL 2 
genes_class<-gur_parsed_df[,c(3,12)] #protein level 2
genes_class <- unique(genes_class)
genes_class<-genes_class[-c(6835,10211,11843),] #OR2A25,GPX6,MKKS remove el no family que es una proteina desconocida con el mismo nomnbre de lo original
rownames(genes_class)<-genes_class[,1]
genes_class <- genes_class[, c("level2_name"), drop = FALSE]

### 1c) PASAR AMBOS NIVELES
both_genes_class<-gur_parsed_df[,c(3,12,14)] #20364
both_genes_class <- unique(both_genes_class) #20345
both_genes_class<-both_genes_class[-c(6835,10211,11843),] #OR2A25,GPX6,MKKS remove el no family que es una proteina desconocida con el mismo nomnbre de lo original
rownames(both_genes_class)<-both_genes_class[,1]
both_genes_class<-both_genes_class[,-1]
genes_class<-both_genes_class #20342
## los que tienen ---- cambiar la palabra por "general"
#both_genes_class<-ifelse(both_genes_class$level2_name)

### parece que solo hemos bajado 2 niveles (no he pasado el nivel 3)

######################### read glow genes matrix ####################

setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")
#2a)original matrix
glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
# voy a probar con la que he rellenado los NAs con el valor máximo + 1
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

#2b)normalized glow matrix
norm_glow_matrix<- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix_normalized.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
glowmatrix<-norm_glow_matrix

###3) assign proteins families to genes
glowmatrix_w_class <- merge(glowmatrix, genes_class, by.x = "row.names", by.y="row.names",all.x = TRUE)
summary(glowmatrix_w_class$level1_name)
unique(glowmatrix_w_class$level1_name)
rownames(glowmatrix_w_class)<-glowmatrix_w_class[,1]
glowmatrix_w_class<-glowmatrix_w_class[,-1]

##mis categorias cubren 20589 genes, la matriz de glow tiene 24757, los genes de glow que no tienen info de protein class se quedan con NA
unique(glowmatrix_w_class$level1_name)



#####
# glowmatrix_w_class <- glowmatrix_w_class[complete.cases(glowmatrix_w_class), ]
# pca_result <- prcomp(glowmatrix_w_class[, -ncol(glowmatrix_w_class)], scale. = TRUE) #perform pca on all columns except last one
# basic<-ggplot2::ggplot(data = as.data.frame(pca_result$x), aes(x = PC1, y = PC2, color = glowmatrix_w_class$level1_name)) +
#   geom_point() +
#   labs(title = "PCA Plot", x = "PC1", y = "PC2")



#4 ) remove no-family genes

#protein level 1
no_fam<-glowmatrix_w_class%>%
  dplyr::filter(level1_name!="NO FAMILY") #nos quedamos con 13557/24757 de genes de glow con info del level 1 de la proteina
unique(no_fam$level1_name)

######## esto no lo uso, solo filtro en base al NO FAMILY de level 1 (protein_level1) #######
#protein level 2: quitar blancos y con ----: me quedo con 15295/24757
# no_fam<-glowmatrix_w_class%>%
#   dplyr::filter(!level2_name %in% c(unique(genes_class$level2_name)[2],unique(genes_class$level2_name)[5]))
# 
# ##creo que habria que quitar los NAs, añadido del 4 septiembre
# no_fam<-no_fam%>%
#   dplyr::filter(!is.na(no_fam$level2_name))
# ####

#### uso esto
#################################### PERFORM TSNE:PINTAR CLASES DE PROTEINAS EN LA PCA #################################### 
# Remove duplicate rows
glowmatrix_w_class<-no_fam #filtrar el no fam por level 1 or level 2. tambien los gene class iniciales meterle las dos categorias de nivel 1 y nivel 2 
#-> al final uso el protein level 1, filtrado por quitar los NOFAMILY, vienen a ser: 13557/24757 genes que tiene glow

tsne_data_original <- glowmatrix_w_class[!duplicated(glowmatrix_w_class), ]
tsne_data_original <- na.omit(tsne_data_original) #13115/13557 quitar los NAs de toda la fila (NA en todos los paneles de glow)

labels_protein_fam_level2<-tsne_data_original$level2_name #13115/13557 
labels_protein_fam_level1<-tsne_data_original$level1_name #13115/13557
# Remove the label column for t-SNE computation
tsne_data <- tsne_data_original[, -ncol(tsne_data_original)]  # Assuming the label is the last column
tsne_data <- tsne_data_original[, -c(ncol(tsne_data_original),ncol(tsne_data_original)-1)]  # quitar las columnas de level 1 y level 2 para el calculo del tsne
# Remove duplicate rows
tsne_data <- tsne_data[!duplicated(tsne_data), ] # esto no hace falta 

# Perform t-SNE
set.seed(42)  # Setting seed for reproducibility
tsne_result <- Rtsne::Rtsne(tsne_data, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300)

# Create a data frame with the t-SNE results
tsne_df <- as.data.frame(tsne_result$Y)

# Add the label back for visualization
# Note: You need to ensure the labels align with the filtered data
# If tsne_data had duplicates removed, you should also filter the labels accordingly -> esto no pasa porque tsne_Data no quita lineas
#level 1
labels_filtered <- tsne_data_original$level1_name[!duplicated(tsne_data_original[, -ncol(tsne_data_original)])]
tsne_df$protein_fam_level1 <- labels_filtered
#level 2
labels_filtered <- tsne_data_original$level2_name[!duplicated(tsne_data_original[, -ncol(tsne_data_original)])]
tsne_df$protein_fam_level2 <- labels_filtered

# Plot t-SNE results using ggplot2 -> solo level 1
library(ggplot2)
##plot todas las familias en el tsne
basic<-ggplot2::ggplot(tsne_df, aes(x = V1, y = V2, color = protein_fam_level1)) +
  geom_point() +
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() +  theme(legend.position = "bottom")   # Adjust legend position as needed
basic
#plot de las familias spliteadas por clase
new <- basic + facet_wrap(~protein_fam_level1, nrow = 7) +  
  geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted")    # Add a horizontal line at y = 0
new

# Extract axis limits from the basic plot
x_limits <- ggplot_build(new)$layout$panel_params[[1]]$x.range
y_limits <- ggplot_build(new)$layout$panel_params[[1]]$y.range

# Filter the dataframe for a specific protein_fam_level1 value -> ejemplo de sacar una categoria del plot en su nivel 1 y pintandola por nivel 2
specific_protein_fam_level1 <- "metabolite interconversion enzyme"  # Replace with the actual value you're interested in
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df %>% filter(protein_fam_level1 == specific_protein_fam_level1)
specific_facet_plot <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level2)) +
  geom_point() +
  labs(title = paste("t-SNE Plot:", specific_protein_fam_level1), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +  # Add a horizontal line at y = 0
  theme_minimal() +
  theme(legend.position = "bottom") +  # Adjust legend position as needed
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits



# more classes:
desired_protein_fam_level1 <-unique(tsne_df$protein_fam_level1)[c(4,8,13,7)]

desired_protein_fam_level1 <-unique(tsne_df$protein_fam_level1)[c(1:6)]
desired_protein_fam_level1 <-unique(tsne_df$protein_fam_level1)[c(7:12)]
desired_protein_fam_level1 <-unique(tsne_df$protein_fam_level1)[c(13:18)]
desired_protein_fam_level1 <-unique(tsne_df$protein_fam_level1)[c(19:24)]
desired_protein_fam_level1 <-unique(tsne_df$protein_fam_level1)[c(25:28)]


plot_list <- lapply(desired_protein_fam_level1, function(level1_value) {
  # Filter the dataframe for the current protein_fam_level1 value
  filtered_df <- tsne_df %>% filter(protein_fam_level1 == level1_value)
  
  # Create the plot for the filtered data with its own legend
  p <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level2)) +
    geom_point() +
    labs(title = paste(level1_value), x = "t-SNE 1", y = "t-SNE 2") +
    geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
    geom_hline(aes(yintercept = 0), linetype = "dotted") +  # Add a horizontal line at y = 0
    theme_minimal() +
    theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
          legend.text = element_text(margin = margin(0, 0, 5, 0)),  # Adjust text margin for better stacking
          legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
          legend.spacing.y = unit(0.2, "cm")) +  # Adjust spacing between legend items
    guides(color = guide_legend(title = NULL))+ # Remove legend title
    coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits
  
  return(p)
})

# Arrange the plots in a grid layout
gridExtra::grid.arrange(grobs = plot_list, ncol = 2)  # Adjust number of columns as needed


###### 









#################################### PERFORM TSNE:PINTAR PANELES EN LA PCA #################################### 
######## aqui lo que queremos ver es como se parecen los paneles en base a sus puntaciones de genes de glow, 
#queremos ver si los paneles del mismo grupo de familias se agrupan y si hay dos paneles de grupos distintos de paneles que se juntan
glowmatrix_wo_class <- glowmatrix[!duplicated(glowmatrix), ]
glowmatrix_wo_class <- na.omit(glowmatrix_wo_class) #19148 -> quitar NAs de glow, info sin clase
tglowmatrix_wo_class<-t(glowmatrix_wo_class)



set.seed(42)  # Setting seed for reproducibility
invert_tsne_result <- Rtsne::Rtsne(tglowmatrix_wo_class, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300)
# Create a data frame with the t-SNE results
invert_tsne_df <- as.data.frame(invert_tsne_result$Y)
invert_tsne_df$panelnames <- rownames(tglowmatrix_wo_class)
invert_tsne_df$panelnames<- gsub("_GA", "", invert_tsne_df$panelnames) #q

### pegar categorias de paneles 
#color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list.txt", header = FALSE, sep = "\t")

colnames(color_panels)<-c("type_panel","panel")
### ajustar paneles
#1) el panel de cancer y el de tumor juntarlos en una categoria
color_panels$type_panel<-ifelse(color_panels$type_panel=="Tumour syndromes","Tumour Syndromes + Cancer Program",color_panels$type_panel)
color_panels$type_panel<-ifelse(color_panels$type_panel=="Cancer Programme","Tumour Syndromes + Cancer Program",color_panels$type_panel)
#2) hay un panel suelto en Haemetological disorders, ponerlo en Haemotoogical and immunological disorders
color_panels$type_panel<-ifelse(color_panels$type_panel=="Haematological and immunological disorders","Haematological disorders",color_panels$type_panel)
color_panels$type_panel<-ifelse(color_panels$type_panel %in% c("Growth disorders",
                                                               "Hearing and ear disorders","Rheumatological disorders",
                                                               "Viral research"), "Other_panels",color_panels$type_panel)


colors_invert_tsne_df <- merge(invert_tsne_df, color_panels, by.x = "panelnames", by.y="panel",all.x = TRUE)
#colors_invert_tsne_df<-colors_invert_tsne_df[-which(colors_invert_tsne_df$panelnames=="Hereditary_neuropathy"),]

# Plot t-SNE results using ggplot2
random_colors_type_panel <- randomcoloR::randomColor(length(unique(colors_invert_tsne_df$type_panel)))

library(ggplot2)
 



#### plot only the points which are similar 
# Load necessary libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# # Define a function to calculate pairwise distances
# calculate_distances <- function(df, threshold) {
#   # Calculate Euclidean distances
#   dist_matrix <- as.matrix(dist(df[, c("V1", "V2")]))
#   
#   # Find pairs of points within the threshold distance
#   close_points <- which(dist_matrix < threshold & dist_matrix > 0, arr.ind = TRUE)
#   
#   # Get unique points to label
#   points_to_label <- unique(c(close_points[,1], close_points[,2]))
#   
#   return(df[points_to_label, ])
# }

# Define a function to calculate pairwise distances and FILTER BY CLASS
calculate_distances <- function(df, threshold) {
  # Calculate Euclidean distances
  dist_matrix <- as.matrix(dist(df[, c("V1", "V2")]))
  
  # Find pairs of points within the threshold distance
  close_points <- which(dist_matrix < threshold & dist_matrix > 0, arr.ind = TRUE)
  
  # Filter pairs by class
  close_points <- close_points[df$type_panel[close_points[,1]] != df$type_panel[close_points[,2]],]
  
  # Get unique points to label
  points_to_label <- unique(c(close_points[,1], close_points[,2]))
  
  return(df[points_to_label, ])
}


# Define the distance threshold (you might need to adjust this)
threshold <- 0.05 # da 5 puntos
threshold <- 0.1 # da 7 puntos
threshold <- 0.2 # da 9 puntos -> esta
threshold <- 0.25 # da 11 puntos
threshold <- 0.3 # da 17 puntos
threshold <- 0.35 # da 20 puntos
# Get points to label
points_to_label <- calculate_distances(colors_invert_tsne_df, threshold)

# Plot with labeled close points
# basic <- ggplot(colors_invert_tsne_df, aes(x = V1, y = V2, color = type_panel)) +
#   geom_point() +
#   ggrepel::geom_text_repel(data = points_to_label, aes(label = panelnames), size = 3, box.padding = 0.35, point.padding = 0.3) +
#   labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2") +
#   theme_minimal()
# basic
# Plot the points, highlighting the ones to label with boxed labels and arrows
#box padding=0.35 original
#point_padding=0.3
basic <- ggplot(colors_invert_tsne_df, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 30) +
  theme_minimal()+ 
  theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2") 
basic



# Generate random colors for each class hue_pal del ggplot2 que es la bonita
set.seed(123)
unique_classes <- unique(colors_invert_tsne_df$type_panel)
num_classes <- length(unique_classes)
random_colors <- sample(hue_pal()(num_classes))
names(random_colors) <- unique_classes
basic<-ggplot2::ggplot(colors_invert_tsne_df, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()+scale_color_manual(values = random_colors)
basic



####################### PREPARAR PLOT PARA EL PAPER #################

#1) el de la pca grande
pca_grande<-ggplot2::ggplot(tsne_df, aes(x = V1, y = V2, color = protein_fam_level1)) +
  geom_point() +
  labs(title = "General Protein families", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal() +   # Adjust legend position as needed
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.2, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.5, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL,nrow = 7, byrow = TRUE)) # Remove legend title
# Extract axis limits from the basic plot
x_limits <- ggplot_build(pca_grande)$layout$panel_params[[1]]$x.range
y_limits <- ggplot_build(pca_grande)$layout$panel_params[[1]]$y.range


#2) los de nivel 1 y 2 del DNA binding transcription factor
specific_protein_fam_level1<-unique(tsne_df$protein_fam_level1)[22] # Replace with the actual value you're interested in
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df %>% filter(protein_fam_level1 == specific_protein_fam_level1)
filtered_df$protein_fam_level2<-ifelse(filtered_df$protein_fam_level2=="--------",
                                       paste0(tools::toTitleCase(specific_protein_fam_level1), " (General)"),
                                       filtered_df$protein_fam_level2)
fam1_level1 <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level1)) +
  geom_point() +
  labs(title = paste(tools::toTitleCase(specific_protein_fam_level1),"family: General category"), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +  # Add a horizontal line at y = 0
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
                           legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
                           legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
                           legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
                           legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL))+ # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits

fam1_level2 <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level2)) +
  geom_point() +
  labs(title = paste(tools::toTitleCase(specific_protein_fam_level1),"family: Specific category"), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL,nrow = 3, byrow = TRUE))+ # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits

#3) los de nivel 2 y 3 del cytoskeletetal protein
specific_protein_fam_level1<-unique(tsne_df$protein_fam_level1)[10]
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df %>% filter(protein_fam_level1 == specific_protein_fam_level1)
filtered_df$protein_fam_level2<-ifelse(filtered_df$protein_fam_level2=="--------",
                                       paste0(tools::toTitleCase(specific_protein_fam_level1), " (General)"),
                                       filtered_df$protein_fam_level2)
fam2_level1 <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level1)) +
  geom_point() +
  labs(title = paste(tools::toTitleCase(specific_protein_fam_level1),"family: General category"), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL)) # Remove legend title

fam2_level2 <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level2)) +
  geom_point() +
  labs(title = paste(tools::toTitleCase(specific_protein_fam_level1),"family: Specific category"), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL,nrow = 4, byrow = TRUE)) + # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits


#3) los de nivel 2 y 3 del protein-modyfiying enzimes
specific_protein_fam_level1<-unique(tsne_df$protein_fam_level1)[7]
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df %>% filter(protein_fam_level1 == specific_protein_fam_level1)
filtered_df$protein_fam_level2<-ifelse(filtered_df$protein_fam_level2=="--------",
                                       paste0(tools::toTitleCase(specific_protein_fam_level1), " (General)"),
                                       filtered_df$protein_fam_level2)
fam3_level1 <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level1)) +
  geom_point() +
  labs(title = paste(tools::toTitleCase(specific_protein_fam_level1),"family: General category"), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL))+ # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits

##se peude poner geom_point(alpha = 0.6)
fam3_level2 <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level2)) +
  geom_point() +
  labs(title = paste(tools::toTitleCase(specific_protein_fam_level1),"family: Specific category"), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL,nrow = 6, byrow = TRUE))+ # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits

#4) paneles genes con etiquetas
plot_paneles <- ggplot(colors_invert_tsne_df, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "vertical",  # Place legend below plot and arrange vertically
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.8,"cm"),
        plot.title = element_text(size = 12)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = NULL,nrow = 5, byrow = TRUE))+ # Remove legend title y poner la leyenda en 3 columnas y 5 filas (15 elementos)
  labs(title = "Gene Panel categories", x = "t-SNE 1", y = "t-SNE 2") 

##### 5) Arrange the plots
library(patchwork)
# Combine fam2, fam1, and fam3 level plots vertically
right_column <- (fam2_level1 + fam2_level2) / 
  (fam1_level1 + fam1_level2) / 
  (fam3_level1 + fam3_level2)

# Combine pca_grande and plot_paneles vertically
left_column <- pca_grande / plot_paneles

# Combine the left and right columns horizontally
combined_plot <- left_column | right_column
combined_plot

ggplot2::ggsave(combined_plot,filename = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/figura1.png",width = 20, height = 11)








################# OTRO FILE DE CLASIFICACION ###########

#https://pantherdb.org/downloads/index.jsp -> PANTHER sequence classification files -> PTHR_18.0_human
PTHR18_human <- read.delim("PTHR18_human.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) #family and genes
colnames(PTHR18_human)<-c("GeneID", "UniProtKB", "GeneName", 
                          "PANTHERFamilySubfamilyID","PANTHERFamilyName",
                          "PANTHERSubfamilyName","MF_GO","BP_GO","CC_GO",
                          "PANTHERClass","PANTHERPathway")
unique_GeneName<-length(unique(PTHR18_human$GeneName)) #19362/19446 genes unicos
unique_PANTHERFamilyName<-length(unique(PTHR18_human$PANTHERFamilyName)) #6971 familias de proteinas
View(as.data.frame(unique(PTHR18_human$PANTHERFamilyName)))

### ARCHIVO CON INFO DE QUÉ HACE CADA FAMILIA 
protein_class <- read.delim("protein_class.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) #family and genes
colnames(protein_class)<-protein_class[4,]
protein_class<-protein_class[-c(1:4),]

#remove empty gene names
non_empty_gene_names <- PTHR18_human %>% 
  dplyr::filter()





