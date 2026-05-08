library(dplyr)
library(tidyr)
library(viridis)
library(ggrepel)
library(patchwork)
library(stringr)


color_panels <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/panels/final_reclasified_panels.tsv", header = FALSE, sep = "\t")
colnames(color_panels)<-c("type_panel","panel")
color_panels$type_panel<-ifelse(color_panels$type_panel %in% c("Growth disorders",
                                                               "Hearing and ear disorders","Rheumatological disorders",
                                                               "Viral research"), "Other_panels",color_panels$type_panel)


glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGenes_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
genes_cluster <- glowmatrix["Cluster"]
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

score <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_score_median.tsv", header = TRUE, sep = ",")


glowmatrix$Symbols <- rownames(glowmatrix)
glowmatrix_w_score <- merge(glowmatrix, score, by = "Symbols")
head(glowmatrix_w_score)

## we plot retinal disorder ranking values (y) for genes with a ranking > 0, vs the median ranking of each gene, coloured by their score value. 

disease <- "Retinal disorders"

#View(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),c("Symbols","Retinal_disorders_GA", "median_ranking", "Score")] %>% dplyr::filter(Retinal_disorders_GA < 250))

basic<-ggplot2::ggplot(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], aes(x = median_ranking, y =Retinal_disorders_GA)) +
  geom_point(aes(color = Score), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5)+
  #scale_color_gradient2(low = "#d68c45", mid = "#cbc9ad", high = "#2c6e49", midpoint = 0.5)+
  #labs(title = "GLOWgenes Ranking", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
  labs(title = "", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
  #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Specificity\nScore")

## we plot retinal disorder ranking values (y), vs the median ranking of each gene, coloured by their score value. 
## this one is not filtering (ranking vlues including 0)
# basic<-ggplot2::ggplot(glowmatrix_w_score, aes(x = median_ranking, y =Retinal_disorders_GA)) +
#   geom_point(aes(color = Score), alpha = 0.7) +
#   scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5)+
#   #scale_color_gradient2(low = "#d68c45", mid = "#cbc9ad", high = "#2c6e49", midpoint = 0.5)+
#   labs(title = "GLOWgenes Ranking", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
#   #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
#   theme_light(base_size = 16) + 
#   #guides(colour=guide_legend(ncol=1)) + 
#   labs(colour = "Score")

basic + geom_text_repel(
  data = subset(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], (Retinal_disorders_GA < 250 & Score > 0.75)),
  aes(label = Symbols),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], (Retinal_disorders_GA < 250 & Score < 0.1)),
  aes(label = Symbols),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()



######3
# Prepare filtered data
filtered_data <- glowmatrix_w_score %>%
  dplyr::filter(Retinal_disorders_GA > 0, Retinal_disorders_GA < 200)

# Top 5 with highest Score
top_high_score <- filtered_data %>%
  arrange(desc(Score)) %>%
  head(10)

# Top 5 with lowest Score
top_low_score <- filtered_data %>%
  arrange(Score) %>%
  head(10)

# Plot
basic + coord_cartesian(clip = "off") +
  geom_label_repel(
    data = top_high_score,
    aes(label = Symbols),
    size = 4, max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    fontface = "bold", fill = "white" , xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)
  ) +
  geom_label_repel(
    data = top_low_score,
    aes(label = Symbols),
    size = 4, max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    fontface = "bold" , fill = "white" , xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)
  ) +
  scale_y_reverse()+ 
  theme(legend.position = c(0.12, 0.3),
        legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))


#####



basic<-ggplot2::ggplot(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], aes(x = Score, y =Retinal_disorders_GA)) +
  geom_point(aes(color = min_ranking)) +
  scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 8000)+
  #scale_color_gradient2(low = "#d68c45", mid = "#cbc9ad", high = "#2c6e49", midpoint = 0.5)+
  #labs(title = "GLOWgenes Ranking", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
  labs(title = "", x = "Score", y = paste("Ranking for ", disease, sep = "" )) +
  #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Min\nRanking")

basic + geom_text_repel(
  data = subset(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], (Retinal_disorders_GA < 250 & Score > 0.75)),
  aes(label = Symbols),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], (Retinal_disorders_GA < 250 & Score < 0.1)),
  aes(label = Symbols),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()




basic<-ggplot2::ggplot(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], aes(x = Score, y =Retinal_disorders_GA)) +
  geom_point(aes(color = median_ranking)) +
  scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 10000)+
  #scale_color_gradient2(low = "#d68c45", mid = "#cbc9ad", high = "#2c6e49", midpoint = 0.5)+
  #labs(title = "GLOWgenes Ranking", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
  labs(title = "", x = "Score", y = paste("Ranking for ", disease, sep = "" )) +
  #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Median\nRanking")

basic + geom_text_repel(
  data = subset(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], (Retinal_disorders_GA < 250 & Score > 0.75)),
  aes(label = Symbols),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], (Retinal_disorders_GA < 250 & Score < 0.1)),
  aes(label = Symbols),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse() #+ ylim(5000,0)



#####################333 Este va a ser el plot bueno

basic<-ggplot2::ggplot(glowmatrix_w_score[which(glowmatrix_w_score$Retinal_disorders_GA>0),], aes(x = Score, y =Retinal_disorders_GA)) +
  geom_point(aes(color = median_ranking)) +
  #scale_color_gradient2(low = "#B95F89", mid = "#CACAAA", high ="#41658A", midpoint = 10000)+
  scale_color_gradient2(low = "#92374D", mid = "#E1BB80", high ="#4A5899", midpoint = 10000)+
  #scale_color_gradient2(low = "#d68c45", mid = "#cbc9ad", high = "#2c6e49", midpoint = 0.5)+
  #labs(title = "GLOWgenes Ranking", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
  labs(title = "", x = "SGDS", y = paste("Ranking for ", disease, sep = "" )) +
  #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Median\nRanking")



######3
# Prepare filtered data
filtered_data <- glowmatrix_w_score %>%
  dplyr::filter(Retinal_disorders_GA > 0, Retinal_disorders_GA < 200)

# Top 5 with highest Score
top_high_score <- filtered_data %>%
  arrange(desc(Score)) %>%
  head(10) 

top_high_score$n_panelapp <- 0
top_high_score[top_high_score$Symbols == "CCDC28B",]$n_panelapp <- 12
top_high_score[top_high_score$Symbols == "TBC1D32",]$n_panelapp <- 12
top_high_score[top_high_score$Symbols == "DTHD1",]$n_panelapp <- 1

# Top 5 with lowest Score
top_low_score <- filtered_data %>%
  arrange(Score) %>%
  head(10)

top_low_score$n_panelapp <- 0
top_low_score[top_low_score$Symbols == "YAP1",]$n_panelapp <- 7
top_low_score[top_low_score$Symbols == "TUBB",]$n_panelapp <- 12
top_low_score[top_low_score$Symbols == "PAX6",]$n_panelapp <- 23
top_low_score[top_low_score$Symbols == "GLI3",]$n_panelapp <- 26
top_low_score[top_low_score$Symbols == "TYR",]$n_panelapp <- 12
top_low_score[top_low_score$Symbols == "TTR",]$n_panelapp <- 18
top_low_score[top_low_score$Symbols == "PEX5",]$n_panelapp <- 20
top_low_score[top_low_score$Symbols == "TUBG1",]$n_panelapp <- 8
top_low_score[top_low_score$Symbols == "NSDHL",]$n_panelapp <- 13
top_low_score[top_low_score$Symbols == "CCT5",]$n_panelapp <- 11

# Plot
basic + coord_cartesian(clip = "off") +
  geom_label_repel(
    data = top_high_score,
    aes(label = Symbols),
    size = 4, max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    fontface = "bold", colour = "black", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)
  ) +
  geom_label_repel(
    data = top_low_score,
    aes(label = Symbols),
    size = 4, max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    fontface = "bold" ,colour = "black", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)
  ) +
  scale_fill_continuous(low = "white", high ="#7b9e87"#, guide = guide_colourbar(theme = theme(legend.direction = "horizontal"))
  ) +
  labs(fill = "PanelApp\nNumber") +
  scale_y_reverse()+
  theme(#legend.position = "none", #legend.position = c(0.10, 0.4),
    legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))



################ voy a guardar los genes green de la v1

## green de la última versión que usamos (la 4)
df_green <- read.table(
  file =  "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/panels/Retinal_disorders_GA.csv" ,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)


df_green <- df_green %>% dplyr::filter(grepl("*Green*", evidence)) 


## los verdes de laprimera versión
df_panelapp_dhr_p <- read.table(
  file =  "/home/yolanda/tblab/yolanda/GLOWgenes/berta/annotation_files/panelapp/retinal_disorders//Retinal_disorders_v1.tsv" ,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)

green <- df_panelapp_dhr_p %>%
  filter(grepl("*Green*", Sources...separated.)) 

# write.table(green["Gene.Symbol"], "/home/yolanda/tblab/yolanda/GLOWgenes/berta/annotation_files/panelapp/retinal_disorders/green_genes_v1.txt", 
#             row.names = FALSE, quote = FALSE)

## el ranking con ToppGene
df_glow_prediction_dhr_v1 <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_DHR_V1_panelapp.csv", header = T, stringsAsFactors = F, check.names=F, sep = ",")
df_glow_prediction_dhr_v1 <- merge(df_glow_prediction_dhr_v1, score, by.x = "GeneSymbol", by.y = "Symbols")

green_after <- df_glow_prediction_dhr_v1 %>% dplyr::filter(GeneSymbol %in% df_green$gene_data.gene_symbol)


library(ggrepel)

basic<-ggplot2::ggplot(df_glow_prediction_dhr_v1, aes(x = median_ranking, y =Rank)) +
  geom_point(aes(color = Score), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5)+
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Specificity\nScore")

basic + geom_text_repel(
  data = subset(df_glow_prediction_dhr_v1[which(df_glow_prediction_dhr_v1$Rank>0),], (Rank < 250 & Score > 0.75)),
  aes(label = GeneSymbol),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(df_glow_prediction_dhr_v1[which(df_glow_prediction_dhr_v1$Rank>0),], (Rank < 250 & Score < 0.1)),
  aes(label = GeneSymbol),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()


green_after <- df_glow_prediction_dhr_v1 %>% dplyr::filter(GeneSymbol %in% df_green$gene_data.gene_symbol)


basic<-ggplot2::ggplot(green_after, aes(x = median_ranking, y =Rank)) +
  geom_point(aes(color = Score), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5)+
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Specificity\nScore")

basic + geom_text_repel(
  data = subset(green_after[which(green_after$Rank>0),], (Rank < 250 & Score > 0.6)),
  aes(label = GeneSymbol),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(green_after[which(green_after$Rank>0),], (Rank < 250 & Score < 0.1)),
  aes(label = GeneSymbol),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()



basic + geom_text_repel(
  data = subset(green_after[which(green_after$Rank>0),], (Score > 0.6)),
  aes(label = GeneSymbol),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(green_after[which(green_after$Rank>0),], (Rank < 250 & Score < 0.1)),
  aes(label = GeneSymbol),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()



basic<-ggplot2::ggplot(green_after, aes(x = Score, y =Rank)) +
  geom_point(aes(color = median_ranking)) +
  #scale_color_gradient2(low = "#B95F89", mid = "#CACAAA", high ="#41658A", midpoint = 10000)+
  scale_color_gradient2(low = "#92374D", mid = "#E1BB80", high ="#4A5899", midpoint = 10000)+
  #scale_color_gradient2(low = "#d68c45", mid = "#cbc9ad", high = "#2c6e49", midpoint = 0.5)+
  #labs(title = "GLOWgenes Ranking", x = "Median Ranking", y = paste("Ranking for ", disease, sep = "" )) +
  #theme_light() +  #theme(legend.position = "right")   # Adjust legend position as needed
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Median\nRanking")


# Plot
basic + coord_cartesian(clip = "off") +
  scale_fill_continuous(low = "white", high ="#7b9e87"#, guide = guide_colourbar(theme = theme(legend.direction = "horizontal"))
  ) +
  geom_text_repel(
    data = subset(green_after[which(green_after$Rank>0),], (Score > 0.5)),
    aes(label = GeneSymbol),
    size = 4, max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    fontface = "bold"
  ) +
  labs(fill = "PanelApp\nNumber") +
  scale_y_reverse()+
  theme(#legend.position = "none", #legend.position = c(0.10, 0.4),
    legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))





