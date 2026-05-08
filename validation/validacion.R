library (openxlsx)
library(dplyr)

retnet = read.xlsx(xlsxFile = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/retnet_2.xlsx", sheet = 1, skipEmptyRows = FALSE)


retnet$Gene <- sub("&#10;.*", "", retnet$`Symbols;OMIM.Numbers`)
retnet$Gene <- sub(",.*", "", retnet$Gene)

# retnet <- retnet %>%
#   filter(!grepl("^MT-", Gene)) # vamos a probar quitando los genes mitocondriales
# 

retnet$Years <- gsub("[^0-9]", " ", retnet$References)  # deja solo números
retnet$Years <- gsub("\\s+", " ", retnet$Years)         # limpia espacios extra
retnet$Years <- trimws(retnet$Years)                    # quita espacios extremos

retnet$Years_list <- regmatches(retnet$References,
                                gregexpr("[0-9]+", retnet$References))

retnet$Year_min <- sapply(retnet$Years_list, function(x) {
  nums <- as.numeric(x)
  years_full <- ifelse(nums <= 26, 2000 + nums, 1900 + nums)
  min(years_full)
})

retnet$Year_max <- sapply(retnet$Years_list, function(x) {
  nums <- as.numeric(x)
  years_full <- ifelse(nums <= 26, 2000 + nums, 1900 + nums)
  max(years_full)
})

retnet$Num_pubs <- lengths(retnet$Years_list)

head(retnet)


df_score = read.table("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", 
                      head=TRUE, sep = ",")


retnet <- merge(retnet, df_score, by.x = "Gene", by.y = "Symbols", all.x = TRUE)
nrow(retnet[is.na(retnet$Score),])
retnet <- retnet[!is.na(retnet$Score), ]
#df_genes_score$Version2_Score <- factor(df_genes_score$Version2_Score, levels = c(5,4,3,2,1,0))

nrow(retnet[is.na(retnet$Score),])
subset(retnet, Year_min > Year_max) # no hay ninguno que tenga los años mal 


plot(retnet$Year_max, retnet$Score,
     xlab = "Año última publicacion",
     ylab = "Score",
     main = "Score vs Año última publicacion",
     pch = 19)

plot(retnet$Year_min, retnet$Score,
     xlab = "Año primera publicacion",
     ylab = "Score",
     main = "Score vs Año primera publicacion",
     pch = 19)


plot(retnet$Num_pubs, retnet$Score,
     xlab = "Número de publicaciones",
     ylab = "Score",
     main = "Score vs Número de publicaciones",
     pch = 19)


library(ggplot2)

ggplot(retnet, aes(x = Year_min, y = Score)) +
  geom_point() +
  labs(title = "Score vs Año primera publicacion",
       x = "Año primera publicacion",
       y = "Score") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

ggplot(retnet, aes(x = Year_max, y = Score)) +
  geom_point() +
  labs(title = "Score vs Año última publicacion",
       x = "Año última publicacion",
       y = "Score") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

ggplot(retnet, aes(x = Num_pubs, y = Score)) +
  geom_point() +
  labs(title = "Score vs Número de publicaciones",
       x = "Año Número de publicaciones",
       y = "Score") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


cor.test(retnet$Year_min, retnet$Score, method = "pearson")


retnet$Year_group <- cut(retnet$Year_min,
                         breaks = c(-Inf, 1990, 2000, 2010, 2020, Inf),
                         labels = c("<1990", "1990-2000", "2000-2010", "2010-2020", ">2020"))

# retnet$Year_group <- cut(retnet$Year_min,
#                          breaks = c(-Inf, 2000, 2015, Inf),
#                          labels = c("<2000", "2000-2015",">2015"))

table(retnet$Year_group)
aggregate(Score ~ Year_group, data = retnet, mean)

ggplot(retnet, aes(x = Year_group, y = Score)) +
  geom_boxplot() +
  theme_minimal()




retnet$N_pub_group <- cut(retnet$Num_pubs,
                         breaks = c(-Inf, 3, Inf),
                         labels = c("<3", ">3"))

table(retnet$N_pub_group)
aggregate(Score ~ N_pub_group, data = retnet, mean)

ggplot(retnet, aes(x = N_pub_group, y = Score)) +
  geom_boxplot() +
  theme_minimal()



modelo <- lm(Score ~ Year_min + Num_pubs, data = retnet)
summary(modelo)



###### HPOs

hpos <- read.delim("~/Downloads/genes_to_phenotype.txt", 
                   sep = "\t", 
                   header = TRUE, 
                   quote = "",       # ignora comillas problemáticas
                   stringsAsFactors = FALSE, 
                   fill = TRUE)      # rellena con NA si faltan columnas


hpos_count <- hpos %>%
  group_by(gene_symbol) %>%
  summarise(n_HPOs = n_distinct(hpo_id)) %>%
  arrange(desc(n_HPOs))

head(hpos_count)

retnet <- merge(retnet, hpos_count, by.x = "Gene", by.y = "gene_symbol")


ggplot(retnet, aes(x = n_HPOs, y = Score)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

cor.test(retnet$n_HPOs, retnet$Score, method = "pearson")


retnet$hpo_group <- cut(retnet$n_HPOs,
                          breaks = c(-Inf, 30, 50, Inf),
                          labels = c("1-30","30-50",">50"),
                          right = TRUE)



as.data.frame(table(retnet$hpo_group))
ggplot(retnet, aes(x = hpo_group, y = Score, fill = hpo_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Número de Hpos asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de hpos") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

aggregate(Score ~ hpo_group, data = retnet, 
          FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))




########## Diseases (jensen)

diseases <- read.table("/home/yolanda/Downloads/human_disease_knowledge_filtered.tsv", sep = "\t")

disease_count <- diseases %>%
  group_by(V2) %>%                           # agrupa por gen
  summarise(n_diseases = n_distinct(V3)) %>% # cuenta enfermedades únicas por DOID
  arrange(desc(n_diseases))                  # opcional: ordenar de más a menos

head(disease_count)

retnet <- merge(retnet, disease_count, by.x = "Gene", by.y = "V2")

ggplot(retnet, aes(x = n_diseases, y = Score)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


ggplot(retnet, aes(x = n_diseases, y = n_HPOs)) +
  geom_point() +
  theme_minimal()

retnet$diseases_group <- cut(retnet$n_diseases,
                             breaks = c(0, 1, 2, Inf),
                             labels = c("1","2",">2"),
                             right = TRUE)



as.data.frame(table(retnet$diseases_group))
ggplot(retnet, aes(x = diseases_group, y = Score, fill = diseases_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Número de diseases asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de diseases") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

aggregate(Score ~ diseases_group, data = retnet, 
          FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))



############# orphanet 

library(readr)

orphanet <- "/home/yolanda/Downloads/gene_disease.csv"
df_orpha <- read_csv("/home/yolanda/Downloads/gene_disease.csv")

df_orpha_count <- df_orpha %>%
  group_by(gene) %>%                           # agrupa por gen
  summarise(n_orpha = n_distinct(orpha_code)) %>% # cuenta enfermedades únicas por DOID
  arrange(desc(n_orpha))                  # opcional: ordenar de más a menos

head(df_orpha_count)

retnet <- merge(retnet, df_orpha_count, by.x = "Gene", by.y = "gene")

ggplot(retnet, aes(x = n_orpha, y = Score)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()



retnet$orpha_group <- cut(retnet$n_orpha,
                           breaks = c(0, 1, 2, 3, Inf),
                           labels = c("1","2","3",">3"),
                           right = TRUE)


retnet$orpha_group <- cut(retnet$n_orpha,
                          breaks = c(0, 3, Inf),
                          labels = c("<3",">3"),
                          right = TRUE)

as.data.frame(table(retnet$orpha_group))
ggplot(retnet, aes(x = orpha_group, y = Score, fill = orpha_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Número de orphas asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de orphas") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

aggregate(Score ~ orpha_group, data = retnet, 
          FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))

###################### vamos a probar con Panelapp 

panelapp_DHR <- list.files(
  path = "/home/yolanda/tblab/yolanda/GLOWgenes/berta/annotation_files/panelapp/retinal_disorders/",
  pattern = "Retinal_disorders_v[0-9]+\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

df_green <- data.frame(Gene.Symbol = character(), 
                       panel = character()
)

df_red <- data.frame(Gene.Symbol = character(), 
                     panel = character()
)

for (p in panelapp_DHR){
  panelapp_version <- basename(p)
  panelapp_version <- tools::file_path_sans_ext(panelapp_version)
  # extraer solo el número
  panelapp_version <- gsub(".*_v([0-9]+)$", "\\1", panelapp_version)
  
  df_panelapp_dhr_p <- read.table(
    file = p,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )

  green <- df_panelapp_dhr_p %>%
    filter(grepl("*Green*", Sources...separated.)) 
  
  green$panel <- panelapp_version
  
  df_green <- rbind(df_green, green[c("Gene.Symbol", "panel")])
  
  red <- df_panelapp_dhr_p %>%
    filter(grepl("*Red*", Sources...separated.)) 
  
  red$panel <- panelapp_version
  
  df_red <- rbind(df_red, red[c("Gene.Symbol", "panel")])
  
}

head(df_green)

df_green$panel <- as.numeric(df_green$panel)

df_first_green <- df_green %>%
  group_by(Gene.Symbol) %>%
  slice_min(panel, with_ties = FALSE) %>%
  ungroup()


head(df_first_green)

df_first_green <- merge(df_first_green, df_score, by.x = "Gene.Symbol", by.y = "Symbols")


ggplot(df_first_green, aes(x = panel, y = Score)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()



df_first_green$panelapp_group <- cut(df_first_green$panel,
                               breaks = c(0,1,2,3, Inf),
                               labels = c("1","2","3", ">3"),
                               right = TRUE)

df_first_green$panelapp_group <- cut(df_first_green$panel,
                               breaks = c(0,3, Inf),
                               labels = c("<3", ">3"),
                               right = TRUE)

as.data.frame(table(df_first_green$panelapp_group))

ggplot(df_first_green, aes(x = panelapp_group, y = Score, fill = panelapp_group)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Número de enfermedades asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de enfermedades") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


df_red$panel <- as.numeric(df_red$panel)

df_first_red <- df_red %>%
  group_by(Gene.Symbol) %>%
  slice_min(panel, with_ties = FALSE) %>%
  ungroup()

df_first_red <- merge(df_first_red, df_score, by.x = "Gene.Symbol", by.y = "Symbols")
head(df_first_red)

ggplot(df_first_red, aes(x = panel, y = Score)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()



df_first_red$panelapp_group <- cut(df_first_red$panel,
                                     breaks = c(0,1,2,3, Inf),
                                     labels = c("1","2","3", ">3"),
                                     right = TRUE)

df_first_red$panelapp_group <- cut(df_first_red$panel,
                                     breaks = c(0,3, Inf),
                                     labels = c("<3", ">3"),
                                     right = TRUE)

as.data.frame(table(df_first_red$panelapp_group))

ggplot(df_first_red, aes(x = panelapp_group, y = Score, fill = panelapp_group)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Número de enfermedades asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de enfermedades") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


df_glow_prediction_dhr <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/results_GA/Retinal_disorders_GA/GLOWgenes_prioritization_Random.txt")


df_first_red <- merge(df_first_red, df_glow_prediction_dhr, by.x = "Gene.Symbol", by.y = "V1")

ggplot(df_first_red, aes(x = panel, y = V3, colour = Score)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


################ voy a guardar los genes green de la v1

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


df_glow_prediction_dhr_v1 <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/berta/annotation_files/panelapp/retinal_disorders/GLOWgenes_prioritization_Random.txt")
df_glow_prediction_dhr_v1 <- merge(df_glow_prediction_dhr_v1, df_score, by.x = "V1", by.y = "Symbols")


ggplot(df_glow_prediction_dhr_v1, aes(x = median_ranking, y = V3 ,colour = Score)) +
  geom_point() +
  theme_minimal() + coord_flip()


library(ggrepel)

basic<-ggplot2::ggplot(df_glow_prediction_dhr_v1, aes(x = median_ranking, y =V3)) +
  geom_point(aes(color = Score), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5)+
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Specificity\nScore")

basic + geom_text_repel(
  data = subset(df_glow_prediction_dhr_v1[which(df_glow_prediction_dhr_v1$V3>0),], (V3 < 250 & Score > 0.75)),
  aes(label = V1),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(df_glow_prediction_dhr_v1[which(df_glow_prediction_dhr_v1$V3>0),], (V3 < 250 & Score < 0.1)),
  aes(label = V1),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()


df_green <- read.table(
  file =  "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/panels/Retinal_disorders_GA.csv" ,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)

df_green <- df_green %>% dplyr::filter(grepl("*Green*", evidence)) 

green_after <- df_glow_prediction_dhr_v1 %>% dplyr::filter(V1 %in% df_green$gene_data.gene_symbol)


basic<-ggplot2::ggplot(green_after, aes(x = median_ranking, y =V3)) +
  geom_point(aes(color = Score), alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0.5)+
  theme_minimal(base_size = 16) + 
  #guides(colour=guide_legend(ncol=1)) + 
  labs(colour = "Specificity\nScore")

basic + geom_text_repel(
  data = subset(green_after[which(green_after$V3>0),], (V3 < 250 & Score > 0.6)),
  aes(label = V1),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(green_after[which(green_after$V3>0),], (V3 < 250 & Score < 0.1)),
  aes(label = V1),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()



basic + geom_text_repel(
  data = subset(green_after[which(green_after$V3>0),], (Score > 0.6)),
  aes(label = V1),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
) + geom_text_repel(
  data = subset(green_after[which(green_after$V3>0),], (V3 < 250 & Score < 0.1)),
  aes(label = V1),
  size = 4, max.overlaps = 25,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"),
  fontface = "bold"
)+ scale_y_reverse()



basic<-ggplot2::ggplot(green_after, aes(x = Score, y =V3)) +
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
    data = subset(green_after[which(green_after$V3>0),], (Score > 0.6)),
    aes(label = V1),
    size = 4, max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    fontface = "bold"
  ) +
  labs(fill = "PanelApp\nNumber") +
  scale_y_reverse()+
  theme(#legend.position = "none", #legend.position = c(0.10, 0.4),
    legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))









#############################################33 AHORA VAMOS A COGER TODOS LOS GENES


###################

df_score <- merge(df_score, disease_count, by.x = "Symbols", by.y = "V2")
df_score <- merge(df_score, hpos_count, by.x = "Symbols", by.y = "gene_symbol")
df_score <- merge(df_score, df_orpha_count, by.x = "Symbols", by.y = "gene")

ggplot(df_score, aes(x = n_diseases, y = Score)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


cor.test(df_score$n_diseases, df_score$Score, method = "pearson")


df_score$Diseases_group <- cut(df_score$n_diseases,
                               breaks = c(0,1, 4,  Inf),
                               labels = c("1","2-4", ">4"),
                               right = TRUE)

as.data.frame(table(df_score$Diseases_group))

ggplot(df_score, aes(x = Diseases_group, y = Score, fill = Diseases_group)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Número de enfermedades asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de enfermedades") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


ggplot(df_score, aes(x = n_HPOs, y = Score)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()

cor.test(df_score$n_HPOs, df_score$Score, method = "pearson")



df_score$HPOs_group <- cut(df_score$n_HPOs,
                           breaks = c(0, 30, 80, Inf),
                           labels = c("0-30","30-80",">80"),
                           right = TRUE)


as.data.frame(table(df_score$HPOs_group))
ggplot(df_score, aes(x = HPOs_group, y = Score, fill = HPOs_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Número de HPOs asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de HPOs") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

aggregate(Score ~ HPOs_group, data = df_score, 
          FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))




ggplot(df_score, aes(x = n_orpha, y = Score)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


cor.test(df_score$n_orpha, df_score$Score, method = "pearson")


df_score$orpha_group <- cut(df_score$n_orpha,
                               breaks = c(0,1, 4,  Inf),
                               labels = c("1","2-4", ">4"),
                               right = TRUE)

as.data.frame(table(df_score$orpha_group))

ggplot(df_score, aes(x = orpha_group, y = Score, fill = orpha_group)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Número de enfermedades asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de enfermedades") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")



##### claro, faltan las predicciones... pruebo con el canal de text mining de diseases

df_score = read.table("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", 
                      head=TRUE, sep = ",")

diseases_textmining <- read.table("/home/yolanda/Downloads/human_disease_textmining_filtered.tsv", sep = "\t")

disease_tm_count <- diseases_textmining %>%
  group_by(V2) %>%                           # agrupa por gen
  summarise(n_diseases_TM = n_distinct(V3)) %>% # cuenta enfermedades únicas por DOID
  arrange(desc(n_diseases_TM))                  # opcional: ordenar de más a menos

df_score <- merge(df_score, disease_tm_count, by.x = "Symbols", by.y = "V2")


ggplot(df_score, aes(x = n_diseases_TM, y = Score)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


cor.test(df_score$n_diseases_TM, df_score$Score, method = "pearson")


df_score$Diseases_TM_group <- cut(df_score$n_diseases_TM,
                               breaks = c(0,1, 4,  Inf),
                               labels = c("1","2-4", ">4"),
                               right = TRUE)

as.data.frame(table(df_score$Diseases_TM_group))

ggplot(df_score, aes(x = Diseases_TM_group, y = Score, fill = Diseases_TM_group)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Número de enfermedades asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de enfermedades") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")




quantiles <- quantile(df_score$n_diseases_TM, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
quantiles

# Crear variable categórica usando cut con los cuartiles
df_score$N_TM_quartiles <- cut(df_score$n_diseases_TM,
                                    breaks = quantiles,
                                    include.lowest = TRUE,
                                    labels = c("Q1","Q2","Q3","Q4"))

as.data.frame(table(df_score$N_TM_quartiles))

ggplot(df_score, aes(x = N_TM_quartiles, y = Score, fill = N_TM_quartiles)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.1) +  # muestra puntos individuales
  labs(x = "Número de Gos",
       y = "SGDS",
       title = "Distribución del Score según número de GOs") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

labels_quartiles <- paste0(
  round(head(quantiles, -1), 2),
  " - ",
  round(tail(quantiles, -1), 2)
)

df_score$N_TM_quartiles <- cut(
  df_score$n_diseases_TM,
  breaks = quantiles,
  include.lowest = TRUE,
  labels = labels_quartiles
)

ggplot(df_score, aes(x = N_TM_quartiles, y = Score, fill = N_TM_quartiles)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, alpha = 0.1) +  # muestra puntos individuales
  labs(x = "N diseases",
       y = "SGDS",
       title = "SGDS distribution depending on N Diseases",
       fill = "N Diseases") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")
















### vamos a quitar a diseases las que estén en la tabla de knowledge

diseases <- read.table("/home/yolanda/Downloads/human_disease_knowledge_filtered.tsv", sep = "\t")

head(diseases_textmining)
head(diseases)

diseases_textmining %>%
  inner_join(diseases, by = c("V2", "V3"))

diseases_textmining_filtrado <- diseases_textmining %>%
  anti_join(diseases, by = c("V2", "V3"))



disease_tm_filtered_count <- diseases_textmining_filtrado %>%
  group_by(V2) %>%                           # agrupa por gen
  summarise(n_diseases_TM_filtered = n_distinct(V3)) %>% # cuenta enfermedades únicas por DOID
  arrange(desc(n_diseases_TM_filtered))                  # opcional: ordenar de más a menos

df_score <- merge(df_score, disease_tm_filtered_count, by.x = "Symbols", by.y = "V2")


ggplot(df_score, aes(x = n_diseases_TM_filtered, y = Score)) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


cor.test(df_score$n_diseases_TM_filtered, df_score$Score, method = "pearson")


df_score$Diseases_TM_filtered_group <- cut(df_score$n_diseases_TM_filtered,
                                  breaks = c(0, 10, 50, Inf),
                                  labels = c("<10","10-50",">50"),
                                  right = TRUE)

as.data.frame(table(df_score$Diseases_TM_filtered_group))

ggplot(df_score, aes(x = Diseases_TM_filtered_group, y = Score, fill = Diseases_TM_filtered_group)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Número de enfermedades asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de enfermedades") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")



# Calcular los cuartiles
quantiles <- quantile(df_score$n_diseases_TM_filtered, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# Crear variable categórica usando cut con los cuartiles
df_score$diseases_TM_quartiles <- cut(df_score$n_diseases_TM_filtered,
                                     breaks = quantiles,
                                     include.lowest = TRUE,
                                     labels = c("Q1","Q2","Q3","Q4"))

as.data.frame(table(df_score$diseases_TM_quartiles))

ggplot(df_score, aes(x = diseases_TM_quartiles, y = Score, fill = diseases_TM_quartiles)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2) +  # muestra puntos individuales
  labs(x = "Número de Diseases asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de Diseases (Text mining)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

# Calcular los cuartiles
quantiles <- quantile(df_score$n_HPOs, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# Crear variable categórica usando cut con los cuartiles
df_score$HPOs_group_quartiles <- cut(df_score$n_HPOs,
                                     breaks = quantiles,
                                     include.lowest = TRUE,
                                     labels = c("Q1","Q2","Q3","Q4"))

as.data.frame(table(df_score$HPOs_group_quartiles))

ggplot(df_score, aes(x = HPOs_group_quartiles, y = Score, fill = HPOs_group_quartiles)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Número de HPOs asociadas (agrupado)",
       y = "Score",
       title = "Distribución del Score según número de HPOs") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")





################# OMIM


df_omim_1998 = read.table("~/Downloads/omim_january_1998_new_gene_symbols.txt", 
                          head=FALSE, sep = ",")

df_omim_1998$date = "1998"

df_omim_2000 = read.table("~/Downloads/omim_april_2000_new_gene_symbols.txt", 
                          head=FALSE, sep = ",")

df_omim_2000$date = "2000"

df_omim_2026 = read.table("~/Downloads/omim_2026_new_gene_symbols.txt", head=FALSE, sep = "\t")

df_omim_2026$date = "2026"


df_omim_all <- rbind(df_omim_1998, df_omim_2000, df_omim_2026)


df_score = read.table("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", 
                      head=TRUE, sep = ",")

df_omim_all <- merge(df_omim_all, df_score, by.x = "V1", by.y = "Symbols")
df_omim_all


ggplot(df_omim_all, aes(x = date, y = Score, fill = date)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Año de inclusión en OMIM",
       y = "Score",
       title = "Distribución del Score según año de inclusión en OMIM") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")




df_omim_2000 = read.table("~/Downloads/omim_2000_new_gene_symbols.txt", 
                          head=FALSE, sep = ",")

df_omim_2000$date = "2000"

df_omim_2005 = read.table("~/Downloads/omim_2005_new_gene_symbols.txt", 
                          head=FALSE, sep = ",")

df_omim_2005$date = "2005"

df_omim_2015 = read.table("~/Downloads/omim_2015_new_gene_symbols.txt", head=FALSE, sep = "\t")

df_omim_2015$date = "2015"


df_omim_2025 = read.table("~/Downloads/omim_2025_new_gene_symbols.txt", head=FALSE, sep = "\t")

df_omim_2025$date = "2025"


df_omim_2026 = read.table("~/Downloads/omim_2026_new_gene_symbols.txt", head=FALSE, sep = "\t")

df_omim_2026$date = "2026"

df_omim_all <- rbind(df_omim_2000, df_omim_2005, df_omim_2015, df_omim_2025, df_omim_2026)


df_score = read.table("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", 
                      head=TRUE, sep = ",")

df_omim_all <- merge(df_omim_all, df_score, by.x = "V1", by.y = "Symbols")
df_omim_all


ggplot(df_omim_all, aes(x = date, y = Score, fill = date)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # muestra puntos individuales
  labs(x = "Año de inclusión en OMIM",
       y = "Score",
       title = "Distribución del Score según año de inclusión en OMIM") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")



aggregate(Score ~ date, data = df_omim_all, 
          FUN = function(x) c(mean=mean(x), sd=sd(x), n=length(x)))










### no sé muy bien por qué hay algunos genes en el green que no están en la lista completa.. Pero tampoco salen en la web de genomic england (https://nhsgms-panelapp.genomicsengland.co.uk/panels/307/v8.0)
## aunque si sale en la de panelapp

panelapp_DHR_GREEN <- "/home/yolanda/tblab/yolanda/GLOWgenes/berta/annotation_files/panelapp/retinal_disorders/Retinal_disorders_v8_green.tsv"
df_panelapp_DHR_GREEN <- read.table(
  file = panelapp_DHR_GREEN,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)

length(df_panelapp_DHR_GREEN$Gene.Symbol)
whole_green <- df_panelapp_dhr_p %>% dplyr::filter(Gene.Symbol %in% df_panelapp_DHR_GREEN$Gene.Symbol)
length(whole_green$Gene.Symbol)

