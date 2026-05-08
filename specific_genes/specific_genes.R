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

disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = TRUE)

### rename CLUSTERS FOR PLOT: 
disease_matrix <- disease_matrix %>%
  mutate(Cluster = case_when(
    Cluster == 4 ~ 1,
    Cluster == 3 ~ 2,
    Cluster == 0 ~ 3,
    Cluster == 2 ~ 4,
    Cluster == 1 ~ 5,
    TRUE ~ Cluster # Default case (if needed)
  ))


n_top_bottom <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/n_top_bottom.csv", header = TRUE, sep = ",")

disease_specific <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_panel_specific_panels.csv", header = TRUE, sep = ",")
colnames(disease_specific) <- c("gene", "n_top", "n_bottom","gene_panel")

class_specific <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_class_specific_class.csv", header = TRUE, sep = ",")
colnames(class_specific) <- c("gene", "class", "top_panels_class")

head(disease_specific)
head(class_specific)

disease_specific_freq <- as.data.frame(table(disease_specific$gene_panel))
colnames(disease_specific_freq) <- c("gene_panel", "count")
disease_specific_freq$gene_panel<- gsub("_GA", "", disease_specific_freq$gene_panel) #q

glowgenes_cluster <- as.data.frame(disease_matrix$Cluster)
rownames(glowgenes_cluster) <- disease_matrix$SYMBOL
colnames(glowgenes_cluster) <- "CLUSTER"

glow_disease_specific <- disease_matrix %>% dplyr::filter(SYMBOL %in% disease_specific$gene)
glow_disease_specific_no_cluster <- glow_disease_specific %>% dplyr::select(-Cluster, -SYMBOL)
row_min <- as.numeric(apply(glow_disease_specific_no_cluster, 1, min))
row_median <- as.numeric(apply(glow_disease_specific_no_cluster, 1, median))

disease_specific$min_value <- row_min
disease_specific$median_value <- row_median
disease_specific$cluster <- glow_disease_specific$Cluster

disease_specific$gene_panel<- gsub("_GA", "", disease_specific$gene_panel) #q

disease_specific <- merge(disease_specific, disease_specific_freq, by = "gene_panel")

#install.packages("ggbeeswarm")
# library(ggbeeswarm)
# 
n_top_bottom_info <- merge(n_top_bottom, disease_specific, by.x = "SYMBOL", by.y = "gene", all.x = TRUE)
n_top_bottom_info <- merge(n_top_bottom_info, glowgenes_cluster, by.x = "SYMBOL", by.y = "row.names")
# 
# 
# 
# n_top_bottom_info %>% dplyr::filter( n_top.x == 1) %>%
# ggplot( aes(x=factor(n_top.x), y=n_bottom.x, colour = as.character(CLUSTER))) + 
#   #geom_point(aes(fill = factor(cluster)), color='#404040', shape=21, size=2)+ 
#   geom_beeswarm() + geom_hline(yintercept=164, linetype="dashed", color = "red")
#   #geom_violin()#+ #scale_y_reverse() + 
#   #coord_flip() +
#   #geom_text_repel(data = filter(n_top_bottom, n_top==1 & n_bottom > 180),aes(label=SYMBOL), size = 4, max.overlaps = 100) 
#   #scale_colour_manual(values=colors_palette)+
#   #scale_fill_manual(values=c("#fdfdd9","#666668","#d387af","#fdb6a1","#9770b1"))+
#   #scale_color_viridis(discrete = TRUE, alpha=0.9, option="A") +
#   #ggtitle("Percentage of top position of each gene") + 
#   #labs(colour = "Cluster") + theme_light(base_size = 16) + guides(colour=guide_legend(ncol=5))+
#   #xlab("Gene Ranking Median") +  ylab("% Top Panel") + theme(legend.position = c(0.75, 0.85),
#   #                                                           legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))
# 



colors_palette <- c("#7EBC89","#C1DBB3","#FAEDCA","#F2C078", "#FE5D26")
# Create the histogram
# hist_plot <- n_top_bottom_info %>%
#   group_by(n_top.x) %>%
#   summarise(count = n()) %>%
#   ggplot(aes(x = as.numeric(n_top.x), y = count)) +
#   geom_col(fill = "blue", alpha = 0.7) +
#   theme_minimal(base_size = 15) +
#   labs(
#     title = "Histogram of n_top.x",
#     x = NULL,  # Remove x-axis label for the top plot
#     y = "Count"
#   ) +
#   theme(axis.text.x = element_blank(),  # Hide x-axis text for alignment
#         axis.ticks.x = element_blank())

n_top_bottom_info_percentage <- n_top_bottom_info %>%
  group_by(n_top.x) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) 

## histogram in %
hist_plot <- n_top_bottom_info %>%
  group_by(n_top.x) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%  # Calculate percentage
  ggplot(aes(x = as.numeric(n_top.x), y = percentage)) +
  geom_col(fill = "#003366", alpha = 0.7) +
  theme_minimal(base_size = 15) +
  labs(
    title = "",
    x = NULL,  # Remove x-axis label for the top plot
    y = "", #"Percentage (%)"
  ) #+ theme(axis.text.x = element_blank(),  # Hide x-axis text for alignment
    #    axis.ticks.x = element_blank())

hist_plot

# Preprocess data: group all values > 20 as "20+"
hist_data <- n_top_bottom_info %>%
  mutate(n_top_group = ifelse(n_top.x > 19, "20+", as.character(n_top.x))) %>%
  group_by(n_top_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    percentage = count / sum(count) * 100,
    n_top_group = factor(n_top_group, levels = c(as.character(0:19), "20+"))  # Keep proper order
  )

# Plot
hist_plot <- ggplot(hist_data, aes(x = n_top_group, y = percentage)) +
  geom_col(fill = "#74a4bc", alpha =1) +
  theme_minimal(base_size = 15) +
  labs(
    title = "",
    x = NULL,
    y = "Percentage (%)"
  )

hist_plot

# Create the jitter plot
jitter_plot <- ggplot(n_top_bottom_info, aes(x = as.numeric(n_top.x), y = n_bottom.x, colour = as.factor(CLUSTER))) +
  geom_jitter(position = position_jitter(0.2), cex = 1.2) +
  scale_colour_manual(values=colors_palette) +
  #geom_hline(yintercept = 164, linetype = "dashed", color = "red") +
  #theme_bw(base_size = 11) +
  theme_minimal(base_size = 15) +
  labs(
    #title = "Jitter Plot of n_bottom.x vs. n_top.x",
    title = "",
    x = "",#"Number of TOP",
    y = "",#"Number of BOTTOM",
    colour = "Cluster"
  ) + guides(colour=guide_legend(ncol=3))+
  theme(legend.position = c(0.70, 0.77), axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), plot.margin = margin(0.1,1,0.5,4, "cm"),
        legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))

jitter_plot
# Combine the two plots
combined_plot <- hist_plot / jitter_plot +
  plot_layout(heights = c(1, 3))  # Adjust relative heights

# Print the combined plot
combined_plot


ggplot(disease_specific, aes(x=reorder(gene_panel, count), y=min_value, color = as.character(cluster), size =n_bottom)) +
  geom_point() +
  theme(panel.background = element_rect(colour="black"), legend.position="right",
        axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
        legend.key.size = unit(0.3, 'cm')) + coord_flip()
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 




ggplot(disease_specific, aes(x=gene_panel, y=min_value, color = median_value, size =n_bottom)) +
  geom_point() + scale_y_reverse() +
  scale_color_gradient(low = "blue", high = "red") +
  theme(panel.background = element_rect(colour="black"), #legend.position="none", legend.position=c(-1,0.9),
        axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
        legend.key.size = unit(0.3, 'cm')) + #coord_flip()
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

ggplot(disease_specific, aes(x=gene_panel, y=min_value, color = n_bottom)) +
  geom_point() + scale_y_reverse() +
  scale_color_gradient(low = "blue", high = "red") +
  theme(panel.background = element_rect(colour="black"), #legend.position="none", legend.position=c(-1,0.9),
        axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
        legend.key.size = unit(0.3, 'cm')) + #coord_flip()
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


ggplot(disease_specific, aes(x=reorder(gene_panel, count), y=min_value, color = n_bottom)) +
  geom_point() + 
  labs(
    #title = "Jitter Plot of n_bottom.x vs. n_top.x",
    title = "",
    y = "Best Ranking Value",
    x = "",
    colour = "    Number of\nbottom rankings"
  ) +
  scale_color_gradient(low = "blue", high = "red") +
  #geom_text_repel(data = filter(disease_specific, min_value < 50 ), aes(label=gene))  +
  theme(panel.background = element_rect(colour="black"), #legend.position="none", legend.position=c(-1,0.9),
        axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
        legend.key.size = unit(0.3, 'cm')) + coord_flip() 

## we add labels

disease_specific$panelnames<- gsub("_", " ", disease_specific$gene_panel) #q

# ggplot(disease_specific, aes(y = reorder(panelnames, count), x = min_value, color = n_bottom)) +
#   geom_point() +
#   scale_color_gradient(low = "blue", high = "red") +
#   geom_text(
#     data = filter(disease_specific, min_value < 50 & n_bottom > 175),
#     aes(label = gene),
#     position = position_jitter(width = 1.2, height = 1.3),
#     hjust = 0,  alpha = 0.7
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 10),
#     axis.title.y = element_text(size = 12)
#   )

## los que su best position es < 20 y los que solo tienen un gen específico
ggplot(disease_specific, aes(y = reorder(panelnames, count), x = min_value, color = n_bottom)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  coord_cartesian(clip = "off") +
  geom_text_repel (
    data = filter(disease_specific, min_value < 50 & n_bottom > 170),
    aes(label = gene), box.padding =  0.5, point.padding = 1, min.segment.length	= 0.1, ylim = c(-Inf, Inf)
  ) +
  labs(
    #title = "Jitter Plot of n_bottom.x vs. n_top.x",
    title = "",
    x = '',#"Best Ranking Value",
    y = "",
    colour = "    Number of\nbottom rankings"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x= element_text(size = 15)
  ) + 
  theme(
    legend.margin = margin(0, 0, 0, 0), # turned off for alignment
    legend.position = c(-0.90, 0.5),
    legend.key.size = unit(1, 'cm'),
    legend.title = element_text(size=12), #change legend title font size
    legend.text = element_text(size=10)
  )


### ponemos el label para solo los genes STBD1,  RIPOR3, STMND1

# ggplot(disease_specific, aes(y = reorder(panelnames, count), x = min_value, color = n_bottom)) +
#   geom_point() +
#   scale_color_gradient(low = "blue", high = "red") +
#   geom_text_repel (
#     data = filter(disease_specific, panelnames %in% c("Corneal dystrophy",  "Corneal abnormalities", "Neurotransmitter disorders", "Glaucoma (developmental)", "Optic neuropathy" )),
#     aes(label = gene)
#   ) +
#   geom_text_repel (
#     data = filter(disease_specific, gene %in% c("STBD1",  "RIPOR3", "STMND1")),
#     aes(label = gene), max.overlaps = 100
#   ) +
#   labs(
#     #title = "Jitter Plot of n_bottom.x vs. n_top.x",
#     title = "",
#     x = "Best Ranking Value",
#     y = "",
#     colour = "    Number of\nbottom rankings"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 10),
#     axis.title.y = element_text(size = 12)
#   )

###################################


class_specific_freq <- as.data.frame(table(class_specific$class))


ggplot(class_specific_freq, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = NULL,
       y = NULL,
       fill = "Gene Type") +
  #scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


ggplot(class_specific_freq, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "#7ebd89") +
  #coord_polar(theta = "y") +
  geom_text(aes(label = Freq), hjust = -0.05, size = 4, color = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal(base_size = 15) +
  labs(title = "Number of Class-specific genes",
       y = "Number of genes",
       x = NULL,
       fill = "Class") + coord_flip()

## in percentage

# Calculate percentage
class_specific_freq <- class_specific_freq %>%
  mutate(Percentage = (Freq / sum(Freq)) * 100)

# Plot with percentages
ggplot(class_specific_freq, aes(x = reorder(Var1, Percentage), y = Percentage)) +
  geom_bar(stat = "identity", fill = "#7ebd89") +
  theme_minimal(base_size = 15) +
  labs(
    title = "Percentage of Class-specific Genes",
    x = NULL,
    y = "Percentage",
    fill = "Class"
  ) +
  coord_flip()



############# Análisis manual

glow_class_specific <- disease_matrix %>% dplyr::filter(SYMBOL %in% class_specific$gene)
glow_class_specific_no_cluster <- glow_class_specific %>% dplyr::select(-Cluster, -SYMBOL)
gene_class_median <- as.numeric(apply(glow_class_specific_no_cluster, 1, median))
gene_class_min <- as.numeric(apply(glow_class_specific_no_cluster, 1, min))

View(cbind(glow_class_specific$SYMBOL,gene_class_min, gene_class_median))
specific_class_min_median <- as.data.frame(cbind(glow_class_specific$SYMBOL,glow_class_specific$Cluster, as.numeric(gene_class_min), gene_class_median, class_specific$class, class_specific$top_panels_class))

View(specific_class_min_median)
View(t(disease_matrix %>% dplyr::filter(SYMBOL == "UQCR10")))
## los genes de lucía
ciliopathies_genes <- read.table("/home/yolanda/Downloads/Panel_Ciliopatias")

class_specific %>% dplyr::filter(gene %in% ciliopathies_genes$V1)
specific_class_min_median %>% dplyr::filter(V1 %in% ciliopathies_genes$V1)

View(t(disease_matrix %>% dplyr::filter(SYMBOL == "TMEM218")))
View(t(disease_matrix %>% dplyr::filter(SYMBOL == "DRC1")))
View(t(disease_matrix %>% dplyr::filter(SYMBOL == "RPUSD1")))

disease_specific %>% dplyr::filter(gene %in% ciliopathies_genes$V1)

#write.table(specific_class_min_median %>% dplyr::filter(V1 %in% ciliopathies_genes$V1), "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_ciliopatias_lucia.tsv", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(specific_class_min_median, "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_especificos_clase_final.tsv", row.names = FALSE, sep = "\t", quote = FALSE)


