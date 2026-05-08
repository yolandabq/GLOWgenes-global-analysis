### Plot de los tsne por panel. Antes y después de Glow

setwd("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/protein_families/")
library(dplyr)
library(factoextra)
library(FactoMineR)
library(scales)
library(randomcoloR)
library(ggplot2)
library(ggrepel)
library(gridExtra)

##########
###

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


#######

### pegar categorias de paneles 
#color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")

colnames(color_panels)<-c("type_panel","panel")

color_panels$type_panel<-ifelse(color_panels$type_panel %in% c("Growth disorders",
                                                               "Hearing and ear disorders","Rheumatological disorders",
                                                               "Viral research"), "Other_panels",color_panels$type_panel)



glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
#glowmatrix_withNA <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
# voy a probar con la que he rellenado los NAs con el valor máximo + 1
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

tglowmatrix_wo_class<-t(glowmatrix)


set.seed(42)  # Setting seed for reproducibility
invert_tsne_result <- Rtsne::Rtsne(tglowmatrix_wo_class, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 300)
# Create a data frame with the t-SNE results
invert_tsne_df <- as.data.frame(invert_tsne_result$Y)
invert_tsne_df$panelnames <- rownames(tglowmatrix_wo_class)
invert_tsne_df$panelnames<- gsub("_GA", "", invert_tsne_df$panelnames) #q


colors_invert_tsne_df_after <- merge(invert_tsne_df, color_panels, by.x = "panelnames", by.y="panel",all.x = TRUE)




# Define the distance threshold (you might need to adjust this)
#threshold <- 0.05 # da 5 puntos
#threshold <- 0.1 # da 7 puntos
threshold <- 0.2 # da 9 puntos -> esta
#threshold <- 0.25 # da 11 puntos
#threshold <- 0.3 # da 17 puntos
#threshold <- 0.35 # da 20 puntos
# Get points to label
points_to_label_after <- calculate_distances(colors_invert_tsne_df_after, threshold)


library("scico")
basic <- ggplot(colors_invert_tsne_df_after, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label_after, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#515153") +
  #geom_label_repel(data = colors_invert_tsne_df_after[colors_invert_tsne_df_after$V1 < -15,], aes(label = panelnames), size = 2.25, box.padding = 1.2, 
  #                 point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  theme_minimal()+ #scale_colo3r_scico_d(palette = "hawaii")   +
  theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot after GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic + xlim(-12, 10) + ylim(-10, 12) #xlim(-32, 10) + ylim(-75, 25)


###


plot_group <- function(group_num) {
  
  colors_invert_tsne_df_after %>%
    mutate(highlight = ifelse(type_panel == group_num, as.character(type_panel), 'Other')) %>%
    ggplot(aes(x = V1, y = V2, color = highlight, alpha = highlight)) +
    geom_point() + 
    scale_color_manual(values = c(setNames("#153A51", group_num), 'Other' = "#E5F4FF")) +  # Highlight one group in red, others in grey
    scale_alpha_manual(values = c(setNames(1, group_num), 'Other' = 0)) +  # Set transparency for "Other"
    labs(title = paste("After GLOW", group_num, "Highlighted"), x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal() + 
    theme(legend.position = "none") + xlim(-32, 10) + ylim(-75, 25) # Remove legend for simplicity
}
# Generate plots for each group and store in a list
plot_list <- list()
for (i in unique(colors_invert_tsne_df_after$type_panel)) {
  plot_list[[i]] <- plot_group(i)
}

# Combine all the plots into a single multiplot using patchwork
combined_plot <- wrap_plots(plot_list) + 
  plot_layout(ncol = 4)  # Adjust the number of columns as per your preference

# Display the combined plot
print(combined_plot)



## PANELES ANTES DE GLOW (PONEMOS UN 1 A TODOS LOS QUE TIENEN UN VALOR Y DEJAMOS 0 A LOS QUE SON SEED)

glowmatrix[glowmatrix>0] <- 1

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


colors_invert_tsne_df_before <- merge(invert_tsne_df, color_panels, by.x = "panelnames", by.y="panel",all.x = TRUE)

points_to_label_before <- calculate_distances(colors_invert_tsne_df_before, threshold)


library("scico")
basic <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label_before, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  #geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$V1 < -15,], aes(label = panelnames), size = 2.25, box.padding = 1.2, 
  #                 point.padding = 0.4,show.legend = FALSE, max.overlaps = 20) +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic + xlim(-12, 10) + ylim(-10, 12)


library("scico")
library(patchwork) 
plot_group <- function(group_num) {
  
  colors_invert_tsne_df_before %>%
    mutate(highlight = ifelse(type_panel == group_num, as.character(type_panel), 'Other')) %>%
    ggplot(aes(x = V1, y = V2, color = highlight, alpha = highlight)) +
    geom_point() + 
    scale_color_manual(values = c(setNames("#153A51", group_num), 'Other' = "#E5F4FF")) +  # Highlight one group in red, others in grey
    scale_alpha_manual(values = c(setNames(1, group_num), 'Other' = 0)) +  # Set transparency for "Other"
    #labs(title = paste("t-SNE plot before GLOW", group_num, "Highlighted"), x = "t-SNE 1", y = "t-SNE 2") +
    labs(title = paste("Before GLOW ", group_num, "Highlighted"), x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal() + 
    theme(legend.position = "none") + xlim(-32, 10) + ylim(-75, 25)  # Remove legend for simplicity
}

# Generate plots for each group and store in a list
plot_list <- list()
for (i in unique(colors_invert_tsne_df_before$type_panel)) {
  plot_list[[i]] <- plot_group(i)
}

# Combine all the plots into a single multiplot using patchwork
combined_plot <- wrap_plots(plot_list) + 
  plot_layout(ncol = 4)  # Adjust the number of columns as per your preference

# Display the combined plot
print(combined_plot)


######################### print both ########

basic_before_before <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label_before, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#101010") +
  #geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$panelnames %in% points_to_label_after$panelnames,], 
  #                 aes(label = panelnames), size = 2.25, box.padding = 1.2, 
  #                 point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + #theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2") + xlim(-12, 10) + ylim(-10, 12)
basic_before_before 

basic_before_after <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  #geom_label_repel(data = points_to_label_before, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
  #                 point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#101010") +
  geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$panelnames %in% points_to_label_after$panelnames,], 
                   aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + #theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2")  + xlim(-12, 10) + ylim(-10, 12)
basic_before_after

#basic_before_before / basic_before_after

basic_after_before <- ggplot(colors_invert_tsne_df_after, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = colors_invert_tsne_df_after[colors_invert_tsne_df_after$panelnames %in% points_to_label_before$panelnames,], 
                   aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#101010") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + 
  labs(title = "t-SNE Plot after GLOW", x = "t-SNE 1", y = "t-SNE 2") + xlim(-12, 10) + ylim(-10, 12)
basic_after_before

#basic_before_before / basic_after_before

basic_after_after <- ggplot(colors_invert_tsne_df_after, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label_after, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + 
  labs(title = "t-SNE Plot after GLOW", x = "t-SNE 1", y = "t-SNE 2") + xlim(-12, 10) + ylim(-10, 12)
basic_after_after 

basic_before_after <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$panelnames %in% points_to_label_after$panelnames,], 
                   aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + #theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2")  + xlim(-12, 10) + ylim(-10, 12)
basic_before_after

grid.arrange(basic_before_before,basic_before_after, 
             basic_after_before,basic_after_after, nrow=2, ncol = 2)


#### sin los límites de eje (para el suplementario)

basic_before_before <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label_before, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#101010") +
  #geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$panelnames %in% points_to_label_after$panelnames,], 
  #                 aes(label = panelnames), size = 2.25, box.padding = 1.2, 
  #                 point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + #theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic_before_before 

basic_before_after <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  #geom_label_repel(data = points_to_label_before, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
  #                 point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#101010") +
  geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$panelnames %in% points_to_label_after$panelnames,], 
                   aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + #theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic_before_after

#basic_before_before / basic_before_after

basic_after_before <- ggplot(colors_invert_tsne_df_after, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = colors_invert_tsne_df_after[colors_invert_tsne_df_after$panelnames %in% points_to_label_before$panelnames,], 
                   aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#101010") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + 
  labs(title = "t-SNE Plot after GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic_after_before

#basic_before_before / basic_after_before

basic_after_after <- ggplot(colors_invert_tsne_df_after, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = points_to_label_after, aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + 
  labs(title = "t-SNE Plot after GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic_after_after 

basic_before_after <- ggplot(colors_invert_tsne_df_before, aes(x = V1, y = V2, color = type_panel)) +
  geom_point() +
  geom_label_repel(data = colors_invert_tsne_df_before[colors_invert_tsne_df_before$panelnames %in% points_to_label_after$panelnames,], 
                   aes(label = panelnames), size = 2.25, box.padding = 1.2, 
                   point.padding = 0.4,show.legend = FALSE, max.overlaps = 20, color = "#878787") +
  theme_minimal()+ #scale_color_scico_d(palette = "hawaii")   +
  theme(legend.position = "none") + #theme(legend.position = "bottom") + 
  labs(title = "t-SNE Plot before GLOW", x = "t-SNE 1", y = "t-SNE 2")
basic_before_after

grid.arrange(basic_before_before,basic_before_after, 
             basic_after_before,basic_after_after, nrow=2, ncol = 2)

















