library(dplyr)
library(factoextra)
library(FactoMineR)
library(scales)
library(randomcoloR)
library(ggplot2)
library(ggrepel)
library(gridExtra)


setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/proteins_families/")
glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)


#1) Read the tab-separated file back into R as a dataframe
gur_parsed_df <- read.table("panther_gene_protein_class_gur.txt", header = TRUE, sep = "\t")
gur_parsed_df <- gur_parsed_df %>%
  dplyr::filter(!GeneSymbol=="unassigned") #20364/20589

### 1c) PASAR AMBOS NIVELES
both_genes_class<-gur_parsed_df[,c(3,12,14)] #20364
both_genes_class <- unique(both_genes_class) #20345
both_genes_class<-both_genes_class[-c(6835,10211,11843),] #OR2A25,GPX6,MKKS remove el no family que es una proteina desconocida con el mismo nomnbre de lo original
rownames(both_genes_class)<-both_genes_class[,1]
both_genes_class<-both_genes_class[,-1]
genes_class<-both_genes_class #20342

######################### read glow genes matrix ####################

#setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")
#2a)original matrix
#glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#glowmatrix_withNA <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
# voy a probar con la que he rellenado los NAs con el valor máximo + 1
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]



###3) assign proteins families to genes
glowmatrix_w_cluster <- merge(glowmatrix, genes_cluster, by = "row.names",all.x = TRUE)
row.names(glowmatrix_w_cluster) <- glowmatrix_w_cluster$Row.names
glowmatrix_w_class <- merge(glowmatrix_w_cluster[-1], genes_class, by = "row.names",all.x = TRUE)

summary(glowmatrix_w_class$level1_name)
unique(glowmatrix_w_class$level1_name)
rownames(glowmatrix_w_class)<-glowmatrix_w_class[,1]
glowmatrix_w_class<-glowmatrix_w_class[,-1]

head(glowmatrix_w_class)

### t-sne by family

tsne_data_all <- glowmatrix_w_cluster[, -ncol(glowmatrix_w_cluster)]  # Assuming the label is the last column
# Perform t-SNE

head(tsne_data_all)
head(unique(tsne_data_all[, -1]))
# al quitar la columna del nombre de los genes, salen algunas filas repetidas. El tsne da problemas así que las quitamos. 

dup_genes <- tsne_data_all[, -1][duplicated(tsne_data_all[, -1]),]
rownames(dup_genes)

tsne_data_all[, -1][rownames(dup_genes),] # son genes que están en un panel de panelapp (Mitochondrial_disorder_with_complex_V_deficiency_GA) y son 0 en el resto. Los quitamos.

set.seed(42)  # Setting seed for reproducibility
#tsne_result <- Rtsne::Rtsne((tsne_data_all[, -1]), dims = 2, perplexity = 50, verbose = TRUE, max_iter = 350, partial_pca=FALSE, check_duplicates = FALSE) # así sale el pez pero boca abajo
tsne_result <- Rtsne::Rtsne((unique(tsne_data_all[, -1])), dims = 2, perplexity = 30, verbose = TRUE, max_iter = 350, partial_pca=FALSE) # así sale en el lado correcto pero no sale muy bien

#duplicated_rows <- rownames(tsne_data_all[duplicated(tsne_data_all[, -1]),])

# Create a data frame with the t-SNE results
tsne_df_all <- as.data.frame(tsne_result$Y)


#write.table(tsne_df_all, file = "tsne_R_glowmatrixfull.tsv", sep = "\t", row.names = FALSE)

#cluster

unique_genes <- row.names(unique(tsne_data_all[, -1]))

tsne_df_all$cluster <- glowmatrix_w_cluster[unique_genes,]$Cluster
tsne_df_all$protein_fam_level1 <- glowmatrix_w_class[unique_genes,]$level1_name
tsne_df_all$protein_fam_level2 <- glowmatrix_w_class[unique_genes,]$level2_name
tsne_df_all$gene <- unique_genes #rownames(glowmatrix_w_class)
head(tsne_df_all)

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


tsne_df_without_nofam <- tsne_df_all %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
# Plot t-SNE results using ggplot2 -> solo level 1
library(ggplot2)
##plot todas las familias en el tsne
basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
  geom_point(aes(colour = protein_fam_level1)) + guides(colour=guide_legend(ncol=1, title = "Protein Family")) +
  labs(title = "t-SNR of protein families", x = "t-SNE 1", y = "t-SNE 2") +
  theme_light(base_size = 16) +  theme(legend.position = "right")   # Adjust legend position as needed
basic

#################################
# a ver, voy a probar a pintarlos por tipo de gen

library(biomaRt)
library(org.Hs.eg.db)  # For human genes; change depending on species
library(GenomicFeatures)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Your gene list (replace with your actual gene names)

gene_list <- tsne_df_all$gene


########## con BiomaRt #####################

# Retrieve information about gene type (coding, non-coding)
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list,
                         mart = ensembl)
# 
# # View the results
# print(gene_annotation)
# 
# nrow(tsne_df_all)
# nrow(merge(tsne_df_all, gene_annotation, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)) ## hay genes que se anotan en dos tipos disintos, por ej: gene_annotation %>% dplyr::filter(hgnc_symbol == "CABIN1")
# tsne_df_all <- merge(tsne_df_all, gene_annotation, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
# 
# tsne_df_without_nofam <- tsne_df_all %>% dplyr::filter(protein_fam_level1 != "NO FAMILY") %>% dplyr::filter(gene_biotype != "protein_coding")
# # Plot t-SNE results using ggplot2 -> solo level 1
# library(ggplot2)
# ##plot todas las familias en el tsne
# basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
#   geom_point(aes(colour = gene_biotype)) + guides(colour=guide_legend(ncol=1, title = "Protein Family")) +
#   labs(title = "t-SNR of protein families", x = "t-SNE 1", y = "t-SNE 2") +
#   theme_light(base_size = 16) +  theme(legend.position = "right")   # Adjust legend position as needed
# basic


################### esto para hacer la pca en vez del tsne
#pca_result <- prcomp(tsne_data_all[, -1], center = TRUE, scale. = TRUE)
#pca_2d <- pca_result$x[, 1:2]  # Take the first two components
#pca_2d <- as.data.frame(pca_2d)
#pca_2d$cluster <- glowmatrix_w_cluster$Cluster
# pca_2d$protein_fam_level1 <- glowmatrix_w_class$level1_name
# pca_2d$protein_fam_level2 <- glowmatrix_w_class$level2_name
# pca_2d$gene <- rownames(glowmatrix_w_class)
# 
# 
# pca_2d$cluster <- factor(pca_2d$cluster, levels = c(4,3,0,2,1))
# 
# basic<-ggplot2::ggplot(pca_2d, aes(x = PC1, y = PC2)) +
#   geom_point(aes(colour = cluster)) + scale_colour_manual(values=colors_palette)+
#   labs(title = "PCA Plot", x = "PC1", y = "PC2") +
#   #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
#   theme_light(base_size = 16) + 
#   theme(legend.position = c(0.90,0.92), legend.background = element_rect(colour = "#797979",linetype="dotdash"),
#         legend.title = element_text(face="bold")) +
#   guides(colour=guide_legend(ncol=5)) + labs(colour = "Cluster")
# basic
# 
# tsne_df_without_nofam <- pca_2d %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")

# Plot t-SNE results using ggplot2 -> solo level 1
# library(ggplot2)
# ##plot todas las familias en el tsne
# basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = PC1, y = PC2)) +
#   geom_point(aes(colour = protein_fam_level1)) +
#   labs(title = "PCA of protein families", x = "PC1", y = "PC2") +
#   theme_light(base_size = 16) +  theme(legend.position = "bottom")   # Adjust legend position as needed
# basic
# #plot de las familias spliteadas por clase
# new <- basic + facet_wrap(~protein_fam_level1, nrow = 7) +  
#   geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
#   geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted")    # Add a horizontal line at y = 0
# new +  theme(legend.position = "none",
#              strip.background=element_rect(colour="black", fill="#ffecce"),
#              strip.text = element_text(colour = 'black'))
# 
# basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = PC1, y = PC2)) +
#   geom_point(aes(colour = protein_fam_level2)) +
#   labs(title = "t-SNE Plot of protein families", x = "t-SNE 1", y = "t-SNE 2") +
#   theme_light(base_size = 16) +  theme(legend.position = "bottom")   # Adjust legend position as needed
# basic
# new <- basic + facet_wrap(~protein_fam_level1, nrow = 7) +  
#   geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
#   geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted")    # Add a horizontal line at y = 0
# new +  theme(legend.position = "none",
#              strip.background=element_rect(colour="black", fill="#ffecce"),
#              strip.text = element_text(colour = 'black'))
# 
# 

# x_limits <- ggplot_build(new)$layout$panel_params[[1]]$x.range
# y_limits <- ggplot_build(new)$layout$panel_params[[1]]$y.range
# # Filter the dataframe for a specific protein_fam_level1 value -> ejemplo de sacar una categoria del plot en su nivel 1 y pintandola por nivel 2
# specific_protein_fam_level1 <- "metabolite interconversion enzyme"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "transporter"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "cytoskeletal protein"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "DNA-binding transcription factor"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "defense/immunity protein"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "gene-specific transcriptional regulator"
# specific_protein_fam_level1 <- "transmembrane signal receptor"
# specific_protein_fam_level1 <- "DNA metabolism protein"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "transfer/carrier protein"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "translational protein"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "intercellular signal molecule"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "protein modifying enzyme"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "RNA metabolism protein"  # Replace with the actual value you're interested in
# 
# 
# # Create the plot for the filtered data with the same axis limits as the original plot
# filtered_df <- tsne_df_without_nofam %>% filter(protein_fam_level1 == specific_protein_fam_level1)
# specific_facet_plot <- ggplot(filtered_df, aes(x = PC1, y = PC2, color = protein_fam_level2)) +
#   geom_point() +
#   labs(title = paste("t-SNE Plot:", specific_protein_fam_level1), x = "t-SNE 1", y = "t-SNE 2") +
#   geom_vline(aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
#   geom_hline(aes(yintercept = 0), linetype = "dotted") +  # Add a horizontal line at y = 0
#   theme_minimal() + labs(colour = "Subfamily") + theme_light(base_size = 16) +
#   theme(legend.position = "bottom") +  # Adjust legend position as needed
#   coord_cartesian(xlim = x_limits, ylim = y_limits)  # Apply the same axis limits
# 
# specific_facet_plot
# 
# specific_facet_plot_split <- specific_facet_plot + facet_wrap(~protein_fam_level2, nrow = 2) +  
#   geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
#   geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted") +   # Add a horizontal line at y = 0
#   theme(legend.position = "none")
# specific_facet_plot_split
# 
# 
# 
#  ################################
# 
# specific_protein_fam_level1 <- "DNA-binding transcription factor"  # Replace with the actual value you're interested in
# specific_protein_fam_level1 <- "DNA metabolism protein"
# specific_protein_fam_level1 <- "transporter"
# specific_protein_fam_level1 <- "defense/immunity protein" 
# specific_protein_fam_level1 <- "translational protein"
# specific_protein_fam_level1 <- "protein modifying enzyme" 
# specific_protein_fam_level1 <- "intercellular signal molecule" 
# specific_protein_fam_level1 <- "RNA metabolism protein"
# specific_protein_fam_level1 <- "cytoskeletal protein"  # Replace with the actual value you're interested in
# 
# filtered_df <- tsne_df_all %>% filter(protein_fam_level1 == specific_protein_fam_level1)
# 
# pca_result <- prcomp(tsne_data_all[filtered_df$gene,][, -1], center = TRUE, scale. = TRUE)
# pca_2d <- pca_result$x[, 1:2]  # Take the first two components
# pca_2d <- as.data.frame(pca_2d)
# pca_2d$cluster <- glowmatrix_w_cluster[filtered_df$gene,]$Cluster
# pca_2d$protein_fam_level1 <- glowmatrix_w_class[filtered_df$gene,]$level1_name
# pca_2d$protein_fam_level2 <- glowmatrix_w_class[filtered_df$gene,]$level2_name
# pca_2d$gene <- rownames(glowmatrix_w_class[filtered_df$gene,])
# 
# pca_2d$cluster <- factor(pca_2d$cluster, levels = c(4,3,0,2,1))
# 
# 
# tsne_df_without_nofam <- pca_2d %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
# # Plot t-SNE results using ggplot2 -> solo level 1
# library(ggplot2)
# ##plot todas las familias en el tsne
# basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = PC1, y = PC2)) +
#   geom_point(aes(colour = protein_fam_level2)) +
#   labs(title = "t-SNE Plot of protein families", x = "t-SNE 1", y = "t-SNE 2") +
#   theme_light(base_size = 16) +  theme(legend.position = "bottom")   # Adjust legend position as needed
# basic
# #plot de las familias spliteadas por clase
# new <- basic + facet_wrap(~protein_fam_level2, nrow = 2) +  
#   geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
#   geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted")    # Add a horizontal line at y = 0
# new +  theme(legend.position = "none",
#              strip.background=element_rect(colour="black", fill="#ffecce"),
#              strip.text = element_text(colour = 'black'))
# 
# 
# 
# family_list <- na.omit(unique(tsne_df_all$protein_fam_level1))
# family_list <- family_list[family_list != "NO FAMILY"]
# family_list <- family_list[family_list != "lyase"]
# family_list <- family_list[family_list != "transferase"]
# family_list <- family_list[family_list != "storage protein"]
# family_list <- family_list[family_list != "viral or transposable element protein"]
# plot_list = list()
# c=1
# for (i in family_list){
#   filtered_df <- tsne_df_all %>% filter(protein_fam_level1 == i)
#   
#   pca_result <- prcomp(tsne_data_all[filtered_df$gene,][, -1], center = TRUE, scale. = TRUE)
#   pca_2d <- pca_result$x[, 1:2]  # Take the first two components
#   pca_2d <- as.data.frame(pca_2d)
#   pca_2d$cluster <- glowmatrix_w_cluster[filtered_df$gene,]$Cluster
#   pca_2d$protein_fam_level1 <- glowmatrix_w_class[filtered_df$gene,]$level1_name
#   pca_2d$protein_fam_level2 <- glowmatrix_w_class[filtered_df$gene,]$level2_name
#   pca_2d$gene <- rownames(glowmatrix_w_class[filtered_df$gene,])
#   
#   pca_2d$cluster <- factor(pca_2d$cluster, levels = c(4,3,0,2,1))
#   
#   
#   tsne_df_without_nofam <- pca_2d %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
#   
#   basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = PC1, y = PC2)) +
#     geom_point(aes(colour = protein_fam_level2)) + 
#     #labs(title = paste0("PCA of protein families: ", i ), x = "t-SNE 1", y = "t-SNE 2") +
#     labs(title = i, x = "", y = "") +
#     theme_light() +  theme(legend.position = "none")
#   
#   #print(basic)
#   plot_list[[c]] = basic
#   c = c+1
# }
# 
# library(patchwork)
# 
# # Combine all plots into a single grid with patchwork
# combined_plot <- wrap_plots(plotlist = plot_list, ncol = 4)  # Adjust ncol to control the grid layout
# 
# # Display the combined plot
# print(combined_plot)

############################################################### hasta aqui con la PCA



##############3 IGUAL PERO CON EL T-SNE #################3


family_list <- na.omit(unique(tsne_df_all$protein_fam_level1))
family_list <- family_list[family_list != "NO FAMILY"]
family_list <- family_list[family_list != "lyase"]
family_list <- family_list[family_list != "transferase"]
family_list <- family_list[family_list != "storage protein"]
family_list <- family_list[family_list != "viral or transposable element protein"]


family_list <- family_list[family_list != "cell junction protein"] # si perplexity es mayor de 20, estas dos no funcionan
family_list <- family_list[family_list != "hydrolase"]


library(patchwork)

plot_list = list()

#for (i in family_list[c(5,7,8,10,15,20)]){
for (i in family_list){  
  filtered_df <- tsne_df_all %>% filter(protein_fam_level1 == i)
  
  set.seed(42)
  tsne_result <- Rtsne::Rtsne(tsne_data_all[filtered_df$gene,][, -1], dims = 2, perplexity = 15, verbose = TRUE, max_iter = 300, partial_pca=FALSE)
  
  tsne_df_family <- as.data.frame(tsne_result$Y)
  tsne_df_family$cluster <- glowmatrix_w_cluster[filtered_df$gene,]$Cluster
  tsne_df_family$protein_fam_level1 <- glowmatrix_w_class[filtered_df$gene,]$level1_name
  tsne_df_family$protein_fam_level2 <- glowmatrix_w_class[filtered_df$gene,]$level2_name
  tsne_df_family$gene <- rownames(glowmatrix_w_class[filtered_df$gene,])
  
  tsne_df_family <- merge(tsne_df_family, gene_annotation, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
  
  tsne_df_family$cluster <- factor(tsne_df_family$cluster, levels = c(4,3,0,2,1))
  
  
  tsne_df_without_nofam <- tsne_df_family %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
  tsne_df_without_nofam[tsne_df_without_nofam$protein_fam_level2 == "--------", "protein_fam_level2"] <- "no subfamily"
  
  basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
    geom_point(data = tsne_df_without_nofam[tsne_df_without_nofam$protein_fam_level2 != "no subfamily",] , aes(colour = protein_fam_level2)) +
    geom_point(data = tsne_df_without_nofam[tsne_df_without_nofam$protein_fam_level2 == "no subfamily",] , colour = "grey", aes(colour = "no subfamily") ) +
    #labs(title = paste0("PCA of protein families: ", i ), x = "t-SNE 1", y = "t-SNE 2") +
    labs(title = i, x = "", y = "") +
    theme_light() +  theme(legend.position = "right") + guides(colour=guide_legend(ncol = 1, title = ""))
  
  # basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
  #   geom_point(aes(colour = gene_biotype)) +
  #   #labs(title = paste0("PCA of protein families: ", i ), x = "t-SNE 1", y = "t-SNE 2") +
  #   labs(title = i, x = "", y = "") +
  #   theme_light() +  theme(legend.position = "bottom") + guides(colour=guide_legend(nrow = 2, title = "")) 
  
  #print(basic)
  plot_list[[i]] = basic
  
}

## sin agrupar
plot_list$`gene-specific transcriptional regulator`
plot_list$`extracellular matrix protein`
plot_list$`metabolite interconversion enzyme`
plot_list$`RNA metabolism protein`

# entre medias
plot_list$`protein modifying enzyme` # esta hay algun grupo que se separa un poco pero esta bastante mezclado
plot_list$`translational protein`
plot_list$`protein-binding activity modulator`
plot_list$`defense/immunity protein`

# agrupados
plot_list$`DNA-binding transcription factor`
plot_list$`intercellular signal molecule`
plot_list$`cytoskeletal protein`
plot_list$transporter



# Combine all plots a single grid with patchwork
combined_plot <- wrap_plots(plotlist = plot_list, nrow = 4)  # Adjust ncol to control the grid layout

# Display the combined plot
print(combined_plot)

# Combine all plots a single grid with patchwork
combined_plot <- wrap_plots(plotlist = plot_list[c(5,7,8,10,15,20)], ncol = 3)  # Adjust ncol to control the grid layout

# Display the combined plot
print(combined_plot)


combined_plot <- wrap_plots(plotlist = plot_list[c("DNA-binding transcription factor", "intercellular signal molecule", 
                                                   "cytoskeletal protein", "transporter", 
                                                   "metabolite interconversion enzyme", "protein modifying enzyme")], nrow = 3)  # Adjust ncol to control the grid layout

combined_plot <- wrap_plots(plotlist = plot_list[c("DNA-binding transcription factor", "intercellular signal molecule")], nrow = 1)  # Adjust ncol to control the grid layout


# Display the combined plot
print(combined_plot)

############3 we add the complete plot
set.seed(42)
tsne_result <- Rtsne::Rtsne(unique(tsne_data_all[, -1]), dims = 2, perplexity = 15, verbose = TRUE, max_iter = 300, partial_pca=FALSE)

tsne_df_family <- as.data.frame(tsne_result$Y)
tsne_df_family$cluster <- glowmatrix_w_cluster[rownames(unique(tsne_data_all[,-1])),]$Cluster
tsne_df_family$protein_fam_level1 <- glowmatrix_w_class[rownames(unique(tsne_data_all[,-1])),]$level1_name
tsne_df_family$protein_fam_level2 <- glowmatrix_w_class[rownames(unique(tsne_data_all[,-1])),]$level2_name
tsne_df_family$gene <- rownames(unique(tsne_data_all[,-1]))

tsne_df_without_nofam <- tsne_df_family %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
tsne_df_without_nofam[tsne_df_without_nofam$protein_fam_level2 == "--------", "protein_fam_level2"] <- "no subfamily"

basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
  geom_point(data = tsne_df_without_nofam[tsne_df_without_nofam$protein_fam_level2 != "no subfamily",] , aes(colour = protein_fam_level1)) +
  #geom_point(data = tsne_df_without_nofam[tsne_df_without_nofam$protein_fam_level2 == "no subfamily",] , colour = "grey", aes(colour = "no subfamily") ) +
  #labs(title = paste0("PCA of protein families: ", i ), x = "t-SNE 1", y = "t-SNE 2") +
  labs(title = i, x = "", y = "") +
  theme_light() +  theme(legend.position = "none") + guides(colour=guide_legend(ncol = 2, title = ""))

basic  ##### esta es la que vamos a coger al final
plot_list[["complete_plot"]] <- basic

#########################33

plot_list[[7]]
new <- plot_list[[7]] + facet_wrap(~protein_fam_level2, nrow = 2) +  
  geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted") +   # Add a horizontal line at y = 0
  theme(legend.position = "none")
new



################################

glowmatrix_w_family <- glowmatrix_w_class[ , !(names(glowmatrix_w_class) %in% c("Cluster","level2_name"))]

median_by_diseases <- glowmatrix_w_family %>% dplyr::group_by(level1_name) %>%
  summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE)))

median_across_diseases <- glowmatrix_w_family %>% 
  rowwise() %>%
  mutate(median_expression = median(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(level1_name) %>%
  summarise(median_expression_across_diseases = median(median_expression, na.rm = TRUE))

write.table(median_across_diseases, file = "median_across_diseases_fam.tsv", sep = "\t", row.names = FALSE)
write.table(median_by_diseases, file = "median_by_diseases_fam.tsv", sep = "\t", row.names = FALSE)

library(tidyr)
library(forcats)

glowmatrix_w_family$gene_name <- row.names(glowmatrix_w_family)
glowmatrix_w_family %>% dplyr::filter(level1_name == "transferase")  %>% dplyr::select(-level1_name) %>%
  pivot_longer(
  cols = -gene_name, # Exclude the gene_family column from reshaping
  names_to = "Disease", # Name of the new column for diseases (x-axis)
  values_to = "Expression" # Name of the new column for values (y-axis)
) %>% dplyr::mutate(Disease = fct_reorder(Disease, Expression, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = Disease, y = Expression)) +
  geom_boxplot() +
  labs(title = "Gene Expression Across Diseases", x = "Disease", y = "Expression") +
  theme_minimal() +
  theme(axis.text.x=element_blank()) # Rotate x-axis labels for readability


data <- glowmatrix_w_family %>% dplyr::filter(level1_name == "transferase")  %>% dplyr::select(-level1_name) %>%
  pivot_longer(
    cols = -gene_name, # Exclude the gene_family column from reshaping
    names_to = "Disease", # Name of the new column for diseases (x-axis)
    values_to = "Expression" # Name of the new column for values (y-axis)
  )

qqnorm(data$Expression)
qqline(data$Expression, col = "red")

# OR Kruskal-Wallis test (non-parametric)
kruskal_result <- kruskal.test(Expression ~ Disease, data = data)
print(kruskal_result)

library("FSA")

dunn_result <- dunnTest(Expression ~ Disease, data = data, method = "bonferroni")
print(dunn_result, dunn.test.results = TRUE)

library(ggpubr)
ggplot(data, aes(x = Disease, y = Expression)) +
  geom_boxplot() +
  stat_compare_means(comparisons = combn(unique(data$Disease), 2, simplify = FALSE),
                     method = "wilcox.test", p.adjust.method = "bonferroni") +
  labs(title = "Pairwise Comparison of Gene Expression Across Diseases",
       x = "Disease", y = "Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


glowmatrix_w_subfamily <- glowmatrix_w_class[ , !(names(glowmatrix_w_class) %in% c("Cluster","level1_name"))]

median_by_diseases <- glowmatrix_w_subfamily %>% dplyr::group_by(level2_name) %>%
  summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE)))

median_across_diseases <- glowmatrix_w_subfamily %>% 
  rowwise() %>%
  mutate(median_expression = median(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(level2_name) %>%
  summarise(median_expression_across_diseases = median(median_expression, na.rm = TRUE))

write.table(median_across_diseases, file = "median_across_diseases_subfam.tsv", sep = "\t", row.names = FALSE)
write.table(median_by_diseases, file = "median_by_diseases_subfam.tsv", sep = "\t", row.names = FALSE)


##########################################

## vamos a hacer los tsne de las familias por separado

specific_protein_fam_level1 <- "DNA-binding transcription factor"  # Replace with the actual value you're interested in
specific_protein_fam_level1 <- "DNA metabolism protein"
specific_protein_fam_level1 <- "transporter"
specific_protein_fam_level1 <- "defense/immunity protein" 
specific_protein_fam_level1 <- "translational protein"
specific_protein_fam_level1 <- "protein modifying enzyme" 
specific_protein_fam_level1 <- "intercellular signal molecule" 

filtered_df <- tsne_df_all %>% filter(protein_fam_level1 == specific_protein_fam_level1)

# Perform t-SNE
set.seed(42)  # Setting seed for reproducibility
tsne_result <- Rtsne::Rtsne(tsne_data_all[filtered_df$gene,][, -1], dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300, partial_pca=FALSE)
tsne_df_family <- as.data.frame(tsne_result$Y)
#write.table(tsne_df_all, file = "tsne_R_glowmatrixfull.tsv", sep = "\t", row.names = FALSE)

#cluster
tsne_df_family$cluster <- glowmatrix_w_cluster[filtered_df$gene,]$Cluster
tsne_df_family$protein_fam_level1 <- glowmatrix_w_class[filtered_df$gene,]$level1_name
tsne_df_family$protein_fam_level2 <- glowmatrix_w_class[filtered_df$gene,]$level2_name
tsne_df_family$gene <- rownames(glowmatrix_w_class[filtered_df$gene,])


tsne_df_without_nofam <- tsne_df_family %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")


basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
  geom_point(aes(colour = protein_fam_level2)) + 
  labs(title = "t-SNE Plot of protein families", x = "t-SNE 1", y = "t-SNE 2") +
  theme_light(base_size = 16) +  theme(legend.position = "bottom")   # Adjust legend position as needed
basic

new <- basic + facet_wrap(~protein_fam_level2, nrow = 2) +  
  geom_vline(data = NULL, aes(xintercept = 0), linetype = "dotted") +  # Add a vertical line at x = 0
  geom_hline(data = NULL, aes(yintercept = 0), linetype = "dotted") +   # Add a horizontal line at y = 0
  theme(legend.position = "none")
new


### hacemos un bucle para plotear todos

family_list <- unique(tsne_df_all$protein_fam_level1)

for (i in family_list){
  filtered_df <- tsne_df_all %>% filter(protein_fam_level1 == i)
  
  # Perform t-SNE
  set.seed(42)  # Setting seed for reproducibility
  tsne_result <- Rtsne::Rtsne(tsne_data_all[filtered_df$gene,][, -1], dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300, partial_pca=FALSE)
  tsne_df_family <- as.data.frame(tsne_result$Y)
  #write.table(tsne_df_all, file = "tsne_R_glowmatrixfull.tsv", sep = "\t", row.names = FALSE)
  
  #cluster
  tsne_df_family$cluster <- glowmatrix_w_cluster[filtered_df$gene,]$Cluster
  tsne_df_family$protein_fam_level1 <- glowmatrix_w_class[filtered_df$gene,]$level1_name
  tsne_df_family$protein_fam_level2 <- glowmatrix_w_class[filtered_df$gene,]$level2_name
  tsne_df_family$gene <- rownames(glowmatrix_w_class[filtered_df$gene,])
  
  
  tsne_df_without_nofam <- tsne_df_family %>% dplyr::filter(protein_fam_level1 != "NO FAMILY")
  
  
  basic<-ggplot2::ggplot(tsne_df_without_nofam, aes(x = V1, y = V2)) +
    geom_point(aes(colour = protein_fam_level2)) + 
    labs(title = "t-SNE Plot of protein families", x = "t-SNE 1", y = "t-SNE 2") +
    theme_light(base_size = 16) +  theme(legend.position = "bottom")
  
  print(basic)
  
}




