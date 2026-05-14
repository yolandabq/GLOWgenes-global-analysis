library(dplyr)
library(tidyr)
library(viridis)
library(ggrepel)
library(patchwork)
library(stringr)
library(gridExtra)
library(cowplot)
library(org.Hs.eg.db)
library(biomaRt)


stat_tops <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", header = TRUE)
disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = TRUE)
rownames(disease_matrix) <- disease_matrix$SYMBOL
tsne <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/tsne.tsv", header = TRUE)
#
gene_types <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/gene_types.csv", header = TRUE, sep = ",")

######
nrow(disease_matrix) # 24757
ncol(disease_matrix) # 211

disease_matrix <- disease_matrix[,2:ncol(disease_matrix)] # we remove the symbol column

disease_matrix %>% dplyr::group_by(Cluster) %>% summarise(n())

mean_genes <- as.data.frame(rowMeans(disease_matrix[, 1:ncol(disease_matrix)-1]))
colnames(mean_genes) <- "Mean"
mean_genes["Cluster"] <- disease_matrix$Cluster

c_order <- mean_genes %>% dplyr::group_by(Cluster) %>%  dplyr::summarise(mean_cluster = mean(Mean)) %>% 
  dplyr::arrange(mean_cluster) %>% dplyr::select(Cluster)


mean_genes["Cluster"] <- factor(disease_matrix$Cluster, levels = c_order$Cluster)
head(mean_genes)

df_mean_genes <- mean_genes
mean_genes["HRAS",]

mean_genes <- mean_genes %>% dplyr::group_by(Cluster) %>% mutate(count_cluster = n())

n_genes_cluster <- as.data.frame(unique(mean_genes[c("Cluster","count_cluster")]))
rownames(n_genes_cluster) <- n_genes_cluster$Cluster
n_genes_cluster$Cluster <- factor(n_genes_cluster$Cluster, levels = c_order$Cluster)

n_genes_cluster <- n_genes_cluster %>% 
  arrange(Cluster) 

mean_genes$count_cluster <- factor(mean_genes$count_cluster, levels =n_genes_cluster$count_cluster)

##### vamos a sacar la mediana: 

# Calculate the median of genes instead of the mean
median_genes <- as.data.frame(apply(disease_matrix[, 1:(ncol(disease_matrix) - 1)], 1, median))
colnames(median_genes) <- "Median"

# Add the Cluster column to the new data frame
median_genes["Cluster"] <- disease_matrix$Cluster

# Group by Cluster, calculate the median for each cluster, and arrange
c_order <- median_genes %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::summarise(median_cluster = mean(Median)) %>% 
  dplyr::arrange(median_cluster) %>% 
  dplyr::select(Cluster)

# Reorder the Cluster factor levels based on the arranged cluster order
median_genes["Cluster"] <- factor(disease_matrix$Cluster, levels = c_order$Cluster)

# Check the first few rows of the result
head(median_genes)

median_genes["HRAS",]

median_mean <- merge(median_genes["Median"], df_mean_genes, by = 'row.names')

median_mean <- median_mean %>% dplyr::mutate(mean_median_diff = Mean - Median)
hist(median_mean$mean_median_diff)

median_mean[order(median_mean$mean_median_diff),]
median_mean[order(-median_mean$mean_median_diff),]

View(median_mean %>% dplyr::filter(mean_median_diff < 10 & mean_median_diff > -10))

# Calculate the median excluding zeros
median_genes <- as.data.frame(
  apply(disease_matrix[, 1:(ncol(disease_matrix) - 1)], 1, function(row) {
    # Remove zeros from the row and compute the median
    median(row[row != 0], na.rm = TRUE)
  })
)
colnames(median_genes) <- "Median"

# Add the Cluster column to the new data frame
median_genes["Cluster"] <- disease_matrix$Cluster

# Group by Cluster, calculate the mean of medians for each cluster, and arrange
c_order <- median_genes %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::summarise(median_cluster = mean(Median, na.rm = TRUE)) %>% 
  dplyr::arrange(median_cluster) %>% 
  dplyr::select(Cluster)

# Reorder the Cluster factor levels based on the arranged cluster order
median_genes["Cluster"] <- factor(disease_matrix$Cluster, levels = c_order$Cluster)

# Check the first few rows of the result
head(median_genes)

median_genes["HRAS",]

