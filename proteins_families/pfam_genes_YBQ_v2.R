########################################### PCA GLOW MATRIX #####################################
setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/proteins_families/")
glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)

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

#setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")
#2a)original matrix
#glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#glowmatrix_withNA <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
# voy a probar con la que he rellenado los NAs con el valor máximo + 1
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")

#2b)normalized glow matrix
#norm_glow_matrix<- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix_normalized.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#glowmatrix<-norm_glow_matrix

###3) assign proteins families to genes
glowmatrix_w_cluster <- merge(glowmatrix, genes_cluster, by = "row.names",all.x = TRUE)
row.names(glowmatrix_w_cluster) <- glowmatrix_w_cluster$Row.names
glowmatrix_w_class <- merge(glowmatrix_w_cluster[-1], genes_class, by = "row.names",all.x = TRUE)

summary(glowmatrix_w_class$level1_name)
unique(glowmatrix_w_class$level1_name)
rownames(glowmatrix_w_class)<-glowmatrix_w_class[,1]
glowmatrix_w_class<-glowmatrix_w_class[,-1]

##mis categorias cubren 20589 genes, la matriz de glow tiene 24757, los genes de glow que no tienen info de protein class se quedan con NA
unique(glowmatrix_w_class$level1_name)


#### t-sne with all genes


tsne_data_all <- glowmatrix_w_cluster[, -c(1,ncol(glowmatrix_w_cluster))]  # Assuming the label is the last column and the first column is the gene name

# Perform t-SNE
set.seed(42)  # Setting seed for reproducibility
tsne_result <- Rtsne::Rtsne(tsne_data_all, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300, partial_pca=TRUE, check_duplicates = FALSE)

# Create a data frame with the t-SNE results
tsne_df_all <- as.data.frame(tsne_result$Y)
#write.table(tsne_df_all, file = "tsne_R_glowmatrixfull.tsv", sep = "\t", row.names = FALSE)


#cluster
tsne_df_all$cluster <- glowmatrix_w_cluster$Cluster
tsne_df_all$protein_fam_level1 <- glowmatrix_w_class$level1_name
tsne_df_all$protein_fam_level2 <- glowmatrix_w_class$level2_name
tsne_df_all$gene <- rownames(glowmatrix_w_class)

colors_palette <- c("#7EBC89","#C1DBB3","#FAEDCA","#F2C078", "#FE5D26")
#colors_palette <- c("#FE5D26","#F2C078","#FAEDCA","#C1DBB3", "#7EBC89")

tsne_df_all$cluster <- factor(tsne_df_all$cluster, levels = c(4,3,0,2,1))
basic<-ggplot2::ggplot(tsne_df_all, aes(x = V1, y = V2)) +
  geom_point(aes(colour = cluster)) + scale_colour_manual(values=colors_palette)+
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2") +
  #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
  theme_light(base_size = 16) + 
  theme(legend.position = c(0.82,0.92), legend.background = element_rect(colour = "#797979",linetype="dotdash"),
        legend.title = element_text(face="bold")) +
  guides(colour=guide_legend(ncol=5)) + labs(colour = "Cluster")
basic


#write.csv(tsne_df_all, "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/tsne_df_allgenes.csv", row.names = FALSE)
#tsne_df_all <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/tsne_df_allgenes.csv", header = TRUE, sep = ",")

tsne_df_without_nofam <- tsne_df_all %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
# Plot t-SNE results using ggplot2 -> solo level 1
library(ggplot2)
##plot todas las familias en el tsne
basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
  geom_point(aes(colour = protein_fam_level1)) +
  labs(title = "t-SNE Plot of protein families", x = "t-SNE 1", y = "t-SNE 2") +
  theme_light(base_size = 16) +  theme(legend.position = "bottom")   # Adjust legend position as needed
basic
#plot de las familias spliteadas por clase
new <- basic + facet_wrap(~protein_fam_level1, nrow = 7) +  
  geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted")    # Add a horizontal line at y = 0
new +  theme(legend.position = "none",
             strip.background=element_rect(colour="black", fill="#ffecce"),
             strip.text = element_text(colour = 'black'))


basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
  geom_point(aes(colour = protein_fam_level2)) + 
  labs(title = "t-SNE Plot of protein families", x = "t-SNE 1", y = "t-SNE 2") +
  theme_light(base_size = 16) +  theme(legend.position = "bottom")   # Adjust legend position as needed
#basic

new <- basic + facet_wrap(~protein_fam_level1, nrow = 7) +  
  geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted") +   # Add a horizontal line at y = 0
  theme(legend.position = "none")
new


x_limits <- ggplot_build(new)$layout$panel_params[[1]]$x.range
y_limits <- ggplot_build(new)$layout$panel_params[[1]]$y.range

# Filter the dataframe for a specific protein_fam_level1 value -> ejemplo de sacar una categoria del plot en su nivel 1 y pintandola por nivel 2
specific_protein_fam_level1 <- "metabolite interconversion enzyme"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "transporter"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "cytoskeletal protein"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "DNA-binding transcription factor"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "defense/immunity protein"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "gene-specific transcriptional regulator"
specific_protein_fam_level1 <- "transmembrane signal receptor"
specific_protein_fam_level1 <- "DNA metabolism protein"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "transfer/carrier protein"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "translational protein"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "intercellular signal molecule"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "protein modifying enzyme"  # Replace with the actual value you're interested in


# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df_without_nofam %>% filter(protein_fam_level1 == specific_protein_fam_level1)
specific_facet_plot <- ggplot(filtered_df, aes(x = V1, y = V2, color = protein_fam_level2)) +
  geom_point() +
  labs(title = paste("t-SNE Plot:", specific_protein_fam_level1), x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +  # Add a horizontal line at y = 0
  theme_minimal() + labs(colour = "Subfamily") + theme_light(base_size = 16) +
  theme(legend.position = "bottom") +  # Adjust legend position as needed
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits

specific_facet_plot

specific_facet_plot_split <- specific_facet_plot + facet_wrap(~protein_fam_level2, nrow = 2) +  
  geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted") +   # Add a horizontal line at y = 0
  theme(legend.position = "none")
specific_facet_plot_split



####################### PREPARAR PLOT PARA EL PAPER #################

#1) el de la pca grande
pca_grande<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2, color = protein_fam_level1)) +
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
specific_protein_fam_level1<-unique(tsne_df_without_nofam$protein_fam_level1)[22] # Replace with the actual value you're interested in
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df_without_nofam %>% filter(protein_fam_level1 == specific_protein_fam_level1)
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
  #labs(title = "", x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_light(base_size = 12) +
  theme(legend.position = "bottom", legend.background = element_rect(colour = "#797979",linetype="dotdash"), 
        legend.box = "vertical",  # Place legend below plot and arrange vertically
        #legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 14)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = "Protein Subfamily",nrow = 4, byrow = TRUE)) + # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits


fam1_level2
#3) los de nivel 2 y 3 del cytoskeletetal protein
specific_protein_fam_level1<-unique(tsne_df_without_nofam$protein_fam_level1)[16]
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df_without_nofam %>% filter(protein_fam_level1 == specific_protein_fam_level1)
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
  #labs(title = "", x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_light(base_size = 12) +
  theme(legend.position = "bottom", legend.background = element_rect(colour = "#797979",linetype="dotdash"), 
        legend.box = "vertical",  # Place legend below plot and arrange vertically
        #legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 14)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = "Protein Subfamily",nrow = 5, byrow = TRUE)) + # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits


fam2_level2


fam1_level2+fam2_level2




#3) los de nivel 2 y 3 del protein-modyfiying enzimes
specific_protein_fam_level1<-unique(tsne_df_without_nofam$protein_fam_level1)[8]
# Create the plot for the filtered data with the same axis limits as the original plot
filtered_df <- tsne_df_without_nofam %>% filter(protein_fam_level1 == specific_protein_fam_level1)
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
  #labs(title = "", x = "t-SNE 1", y = "t-SNE 2") +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme_light(base_size = 12) +
  theme(legend.position = "bottom", legend.background = element_rect(colour = "#797979",linetype="dotdash"), 
        legend.box = "vertical",  # Place legend below plot and arrange vertically
        #legend.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove legend margins
        legend.text = element_text(margin = margin(0, 0, 0, 0)),  # Adjust text margin for better stacking
        legend.key.size = unit(0.5, "lines"),  # Reduce size of legend keys
        legend.spacing.y = unit(0.2, "cm"),
        plot.title = element_text(size = 14)) +  # Adjust spacing between legend items
  guides(color = guide_legend(title = "Protein Subfamily",nrow = 4, byrow = TRUE)) + # Remove legend title
  coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits
fam3_level2

fam1_level2+fam2_level2+fam3_level2










###################33 vemos si hay más o menos enfermedades por protein family 


diseases_lower100 <- apply(glowmatrix[, -1], 1, function(x) sum(x <= 500))
diseases_lower100 <- as.data.frame(diseases_lower100)

diseases_lower100 <- merge(diseases_lower100, genes_class, by="row.names")

diseases_lower100_level1 <- diseases_lower100 %>% dplyr::group_by(level1_name) %>% dplyr::summarise(
  median_lower100 = mean(diseases_lower100)
)
View(diseases_lower100_level1)



##################









diseases_lower100_level2 <- diseases_lower100 %>% dplyr::group_by(level2_name) %>% dplyr::mutate(
  median_lower100_level2 = median(diseases_lower100)
)
View(unique(diseases_lower100_level2[c("level1_name","level2_name","median_lower100_level2")]))

n_genes_fam <- tsne_df_without_nofam %>% dplyr::group_by(level1_name) %>% dplyr::count()

zeros_per_gene <- as.data.frame(rowSums(glowmatrix == 0))
colnames(zeros_per_gene) <- "N_panels_input"
zeros_per_gene$gene <- rownames(zeros_per_gene)

cols_with_zeros <- apply(glowmatrix, 1, function(row) {
  gsub("_GA", "", names(glowmatrix)[which(row == 0)])
})

cols_with_zeros

cols_with_zeros <- lapply(X = cols_with_zeros, FUN = function(t) gsub(pattern = "_GA", replacement = "", x = t, fixed = TRUE))
cols_with_zeros <- lapply(X = cols_with_zeros, FUN = function(t) paste(x = t, collapse = ";;"))
cols_with_zeros_df <- data.frame(
  gene = names(cols_with_zeros),
  panels_input = unname(unlist(cols_with_zeros)),
  stringsAsFactors = FALSE
)

gene_fam_npanels <- merge(tsne_df_without_nofam, zeros_per_gene, by = "gene")
gene_fam_npanels <- merge(gene_fam_npanels, cols_with_zeros_df, by = "gene")

panels_class <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
colnames(panels_class) <- c("Class", "Panel")

panels_class$Class<-ifelse(panels_class$Class %in% c("Growth disorders",
                                                               "Hearing and ear disorders","Rheumatological disorders",
                                                               "Viral research"), "Other_panels",panels_class$Class)


# Create a lookup table from df2
lookup_table <- setNames(panels_class$Class, panels_class$Panel)

# Function to map disorders to categories
map_disorders_to_categories <- function(panels, lookup_table) {
  # Split the disorders by ';;'
  panel_list <- strsplit(panels, ";;")[[1]]
  
  # Look up the categories for each disorder and remove any missing ones
  categories <- unique(lookup_table[panel_list])
  categories <- categories[!is.na(categories)]
  
  # Return the categories as a string separated by ';;'
  return(paste(categories, collapse = ";;"))
}

# Apply the function to the last column of df1
gene_fam_npanels$Mapped_Categories <- sapply(gene_fam_npanels$panels_input, map_disorders_to_categories, lookup_table = lookup_table)

# View the result
gene_fam_npanels

gene_count_df <- gene_fam_npanels %>% dplyr::group_by(level1_name) %>%
  summarise(gene_count_fam = n(), .groups = 'drop')

gene_fam_npanels <- gene_fam_npanels %>%
  left_join(gene_count_df, by = "level1_name")

gene_fam_npanels[(gene_fam_npanels$Mapped_Categories == ""), "Mapped_Categories"] <- NA

unique_cat_df <- gene_fam_npanels %>% dplyr::group_by(level1_name) %>% 
  summarise(
    categories = paste(na.omit(Mapped_Categories), collapse = ";;"),
    unique_categories = paste(unique(strsplit(categories, ";;")[[1]]), collapse = ";;"), 
    n_unique_categories = length(strsplit(unique_categories, ";;")[[1]])
  )

unique_cat_df

gene_fam_npanels <- gene_fam_npanels %>%
  left_join(unique_cat_df[-c(2)], by = "level1_name")


unique(gene_fam_npanels[c("level1_name", "gene_count_fam", "n_unique_categories")])

############### Ahora lo vamos a hacer por nivel 2, a ver si sale algo más claro



n_genes_fam_level2 <- tsne_df_without_nofam %>% dplyr::group_by(level2_name) %>% dplyr::count()

zeros_per_gene <- as.data.frame(rowSums(glowmatrix == 0))
colnames(zeros_per_gene) <- "N_panels_input"
zeros_per_gene$gene <- rownames(zeros_per_gene)

cols_with_zeros <- apply(glowmatrix, 1, function(row) {
  gsub("_GA", "", names(glowmatrix)[which(row == 0)])
})

cols_with_zeros

cols_with_zeros <- lapply(X = cols_with_zeros, FUN = function(t) gsub(pattern = "_GA", replacement = "", x = t, fixed = TRUE))
cols_with_zeros <- lapply(X = cols_with_zeros, FUN = function(t) paste(x = t, collapse = ";;"))
cols_with_zeros_df <- data.frame(
  gene = names(cols_with_zeros),
  panels_input = unname(unlist(cols_with_zeros)),
  stringsAsFactors = FALSE
)

gene_fam_npanels_level2 <- merge(tsne_df_without_nofam, zeros_per_gene, by = "gene")
gene_fam_npanels_level2 <- merge(gene_fam_npanels_level2, cols_with_zeros_df, by = "gene")

panels_class <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
colnames(panels_class) <- c("Class", "Panel")

panels_class$Class<-ifelse(panels_class$Class %in% c("Growth disorders",
                                                     "Hearing and ear disorders","Rheumatological disorders",
                                                     "Viral research"), "Other_panels",panels_class$Class)


# Create a lookup table from df2
lookup_table <- setNames(panels_class$Class, panels_class$Panel)

# Function to map disorders to categories
map_disorders_to_categories <- function(panels, lookup_table) {
  # Split the disorders by ';;'
  panel_list <- strsplit(panels, ";;")[[1]]
  
  # Look up the categories for each disorder and remove any missing ones
  categories <- unique(lookup_table[panel_list])
  categories <- categories[!is.na(categories)]
  
  # Return the categories as a string separated by ';;'
  return(paste(categories, collapse = ";;"))
}

# Apply the function to the last column of df1
gene_fam_npanels_level2$Mapped_Categories <- sapply(gene_fam_npanels_level2$panels_input, map_disorders_to_categories, lookup_table = lookup_table)

# View the result
gene_fam_npanels_level2

gene_count_level2_df <- gene_fam_npanels_level2 %>% dplyr::group_by(level2_name) %>%
  summarise(gene_count_fam = n(), .groups = 'drop')

gene_fam_npanels_level2 <- gene_fam_npanels_level2 %>%
  left_join(gene_count_level2_df, by = "level2_name")


gene_fam_npanels_level2[(gene_fam_npanels_level2$Mapped_Categories == ""), "Mapped_Categories"] <- NA

unique_cat_level2_df <- gene_fam_npanels_level2 %>% dplyr::group_by(level2_name) %>% 
  summarise(
    categories = paste(na.omit(Mapped_Categories), collapse = ";;"),
    unique_categories = paste(unique(strsplit(categories, ";;")[[1]]), collapse = ";;"), 
    n_unique_categories = length(strsplit(unique_categories, ";;")[[1]])
  )

unique_cat_level2_df

gene_fam_npanels_level2 <- gene_fam_npanels_level2 %>%
  left_join(unique_cat_level2_df[-c(2)], by = "level2_name")


n_inputs_level2 <- gene_fam_npanels_level2 %>% dplyr::group_by(level2_name) %>% summarise(
  n_inputs_level2 = sum(N_panels_input)
)

gene_fam_npanels_level2 <- gene_fam_npanels_level2 %>%
  left_join(n_inputs_level2, by = "level2_name")


gene_fam_npanels_level2[(gene_fam_npanels_level2$panels_input == ""), "panels_input"] <- NA

unique_cat_level2_df <- gene_fam_npanels_level2 %>% dplyr::group_by(level2_name) %>% 
  summarise(
    panels_input_level2 = paste(na.omit(panels_input), collapse = ";;"),
    n_panels_input_level2all = length(strsplit(panels_input_level2, ";;")[[1]]),
    panels_input_unique = paste(unique(strsplit(panels_input_level2, ";;")[[1]]), collapse = ";;"), 
    n_panels_input_unique = length(strsplit(panels_input_unique, ";;")[[1]])
  )

gene_fam_npanels_level2 <- gene_fam_npanels_level2 %>%
  left_join(unique_cat_level2_df, by = "level2_name")


View(unique(gene_fam_npanels_level2[c("level1_name" ,"level2_name", "gene_count_fam","n_panels_input_level2all" ,"n_panels_input_unique", "n_unique_categories", "unique_categories")]))

#############33 JACCARD SIMILARITY###########3 en principio nada de esto lo hago ##########

zeros_per_gene <- as.data.frame(rowSums(glowmatrix == 0))
colnames(zeros_per_gene) <- "N_panels_input"
zeros_per_gene$gene <- rownames(zeros_per_gene)

cols_with_zeros <- apply(glowmatrix, 1, function(row) {
  gsub("_GA", "", names(glowmatrix)[which(row == 0)])
})

cols_with_zeros

cols_with_zeros <- lapply(X = cols_with_zeros, FUN = function(t) gsub(pattern = "_GA", replacement = "", x = t, fixed = TRUE))

# Initialize an empty data frame
df <- data.frame(Gene = character(), Panel = character(), stringsAsFactors = FALSE)

# Loop through each gene in the list
for (gene in names(cols_with_zeros)) {
  if (length(cols_with_zeros[[gene]]) == 0) {
    # If the list is empty, add "NO PANEL"
    df <- rbind(df, data.frame(gene = gene, Panel = "NO PANEL"))
  } else {
    # Otherwise, add each panel
    for (panel in cols_with_zeros[[gene]]) {
      df <- rbind(df, data.frame(gene = gene, Panel = panel))
    }
  }
}

colnames(df) <- c("gene", "panel")

merge(df, tsne_df_without_nofam, by=gene)

gene_fam <- tsne_df_without_nofam[c("protein_fam_level1", "protein_fam_level2", "gene")]

gene_fam_pbam <- gene_fam[gene_fam$protein_fam_level1 == "protein-binding activity modulator",]

gene_fam_pbam <- merge(gene_fam_pbam, df, by="gene", all.x = TRUE, all.y = FALSE)



######################3 Plot de los PANELES

#glowmatrix_wo_class <- glowmatrix[!duplicated(glowmatrix), ]
#glowmatrix_wo_class <- na.omit(glowmatrix_wo_class) #19148 -> quitar NAs de glow, info sin clase
#tglowmatrix_wo_class<-t(glowmatrix_wo_class)
tglowmatrix_wo_class<-t(glowmatrix)

#glowmatrix_wo_class <- glowmatrix_withNA[!duplicated(glowmatrix_withNA), ]
#glowmatrix_wo_class <- na.omit(glowmatrix_wo_class) #19148 -> quitar NAs de glow, info sin clase
#tglowmatrix_wo_class<-t(glowmatrix_wo_class)

set.seed(42)  # Setting seed for reproducibility
invert_tsne_result <- Rtsne::Rtsne(tglowmatrix_wo_class, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300)
# Create a data frame with the t-SNE results
invert_tsne_df <- as.data.frame(invert_tsne_result$Y)
invert_tsne_df$panelnames <- rownames(tglowmatrix_wo_class)
invert_tsne_df$panelnames<- gsub("_GA", "", invert_tsne_df$panelnames) #q

### pegar categorias de paneles 
#color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")

colnames(color_panels)<-c("type_panel","panel")
### ajustar paneles
#1) el panel de cancer y el de tumor juntarlos en una categoria
#color_panels$type_panel<-ifelse(color_panels$type_panel=="Tumour syndromes","Tumour Syndromes + Cancer Program",color_panels$type_panel)
#color_panels$type_panel<-ifelse(color_panels$type_panel=="Cancer Programme","Tumour Syndromes + Cancer Program",color_panels$type_panel)
#2) hay un panel suelto en Haemetological disorders, ponerlo en Haemotoogical and immunological disorders
#color_panels$type_panel<-ifelse(color_panels$type_panel=="Haematological and immunological disorders","Haematological disorders",color_panels$type_panel)
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
#install.packages("scico")
library("scico")
basic <- ggplot(colors_invert_tsne_df, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  geom_label_repel(data = colors_invert_tsne_df[colors_invert_tsne_df$V1 < -15,], aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  #theme(legend.position = "bottom") + 
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect(fill = "#d9d9d9", colour = "#797979"))+
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2")+
  guides(color = guide_legend(title = "Panel Class", ncol=2))
basic

colors_invert_tsne_df_nooutlayer <- colors_invert_tsne_df %>% dplyr::filter(panelnames != "Rare_syndromic_craniosynostosis_or_isolated_multisuture_synostosis")
basic <- ggplot(colors_invert_tsne_df_nooutlayer, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  geom_label_repel(data = colors_invert_tsne_df_nooutlayer[colors_invert_tsne_df_nooutlayer$V1 < -15,], aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  #theme(legend.position = "bottom") + 
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_rect(fill = "#d9d9d9", colour = "#797979"))+
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2")+
  guides(color = guide_legend(title = "Panel Class", ncol=2))
basic



setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")
library(patchwork) 


# Function to plot for each group with one highlighted and others in grey
plot_group <- function(group_num) {
  
  colors_invert_tsne_df_nooutlayer %>%
    mutate(highlight = ifelse(type_panel == group_num, as.character(type_panel), 'Other')) %>%
    ggplot(aes(x = V1, y = V2, color = highlight, alpha = highlight)) +
    geom_point() +
    scale_color_manual(values = c(setNames("#153A51", group_num), 'Other' = "#E5F4FF")) +  # Highlight one group in red, others in grey
    scale_alpha_manual(values = c(setNames(1, group_num), 'Other' = 0)) +  # Set transparency for "Other"
    labs(title = paste("t-SNE plot", group_num, "Highlighted"), x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal() + 
    theme(legend.position = "none")  # Remove legend for simplicity
}

# Generate plots for each group and store in a list
plot_list <- list()
for (i in unique(colors_invert_tsne_df_nooutlayer$type_panel)) {
  plot_list[[i]] <- plot_group(i)
}

# Combine all the plots into a single multiplot using patchwork
combined_plot <- wrap_plots(plot_list) + 
  plot_layout(ncol = 4)  # Adjust the number of columns as per your preference

# Display the combined plot
print(combined_plot)




# Function to plot for each group by highlighting type_panel and labeling panel_name
plot_group <- function(group_num) {
  
  colors_invert_tsne_df_nooutlayer %>%
    mutate(highlight = ifelse(type_panel == group_num, as.character(type_panel), 'Other')) %>%
    ggplot(aes(x = V1, y = V2, color = highlight)) +
    geom_point() +
    geom_text_repel(aes(label = ifelse(type_panel == group_num, panelnames, "")), 
                    vjust = -1, hjust = 0.5, size = 3) +  # Add labels for highlighted group using panel_name
    scale_color_manual(values = c(setNames("#2196f3", group_num), 'Other' = "#f0d7db")) +  # Highlight one group in red, others in grey
    labs(title = paste("Group", group_num, "Highlighted")) +
    theme_minimal() +
    theme(legend.position = "none")  # Remove legend for simplicity
}

# Generate plots for each group and save
for (i in unique(colors_invert_tsne_df_nooutlayer$type_panel)) {
  p <- plot_group(i)
  print(p)  # Print the plot
  ggsave(filename = paste0("tsne_group_", i, ".png"), plot = p)  # Save the plot
}


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


###############

library(FactoMineR)

pca2.nci <- PCA(X = tglowmatrix_wo_class, scale.unit = TRUE, ncp = 64, graph = FALSE)


library(factoextra)

fviz_pca_ind(pca2.nci, geom.ind = "point", 
             col.ind = "#FC4E07", 
             axes = c(1, 2), 
             pointsize = 1.5) 

##########









