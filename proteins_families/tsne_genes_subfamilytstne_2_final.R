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

######################### filter glow genes matrix ####################


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

## extraemos las familias 
family_list <- na.omit(unique(glowmatrix_w_class$level1_name))
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
  filtered_df <- rownames(glowmatrix_w_class %>% filter(level1_name == i))
  
  set.seed(42)
  tsne_result <- Rtsne::Rtsne(tsne_data_all[filtered_df,][, -1], dims = 2, perplexity = 15, verbose = TRUE, max_iter = 300, partial_pca=FALSE)
  
  tsne_df_family <- as.data.frame(tsne_result$Y)
  tsne_df_family$cluster <- glowmatrix_w_cluster[filtered_df,]$Cluster
  tsne_df_family$protein_fam_level1 <- glowmatrix_w_class[filtered_df,]$level1_name
  tsne_df_family$protein_fam_level2 <- glowmatrix_w_class[filtered_df,]$level2_name
  tsne_df_family$gene <- rownames(glowmatrix_w_class[filtered_df,])
  
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


















