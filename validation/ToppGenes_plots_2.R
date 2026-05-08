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


stat_tops <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_stat_tops_vfinal.tsv", header = TRUE)
stat_tops_GLOW <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", header = TRUE)
disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGenes_clustered_no_NA.tsv", header = TRUE)
rownames(disease_matrix) <- disease_matrix$SYMBOL
#
#gene_types <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/gene_types.csv", header = TRUE, sep = ",")

### rename CLUSTERS FOR PLOT: 
disease_matrix <- disease_matrix %>%
  mutate(Cluster = case_when(
    Cluster == 1 ~ 1,
    Cluster == 3 ~ 2,
    Cluster == 0 ~ 3,
    Cluster == 2 ~ 4,
    Cluster == 4 ~ 5,
    TRUE ~ Cluster # Default case (if needed)
  ))

stat_tops  <- stat_tops %>%
  mutate(cluster = case_when(
    cluster == 1 ~ 1,
    cluster == 3 ~ 2,
    cluster == 0 ~ 3,
    cluster == 2 ~ 4,
    cluster == 4 ~ 5,
    TRUE ~ cluster # Default case (if needed)
  ))


stat_tops_GLOW  <- stat_tops_GLOW %>%
  mutate(cluster = case_when(
    cluster == 4 ~ 1,
    cluster == 3 ~ 2,
    cluster == 0 ~ 3,
    cluster == 2 ~ 4,
    cluster == 1 ~ 5,
    TRUE ~ cluster # Default case (if needed)
  ))

# gene_types  <- gene_types %>%
#   mutate(cluster = case_when(
#     cluster == 4 ~ 1,
#     cluster == 3 ~ 2,
#     cluster == 0 ~ 3,
#     cluster == 2 ~ 4,
#     cluster == 1 ~ 5,
#     TRUE ~ cluster # Default case (if needed)
#   ))
# 
# 
# write.table(disease_matrix, "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_RENAMED_clustered_no_NA.tsv", row.names = FALSE, quote = FALSE, sep = "\t")


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

# mean_genes %>%
#   ggplot( aes(x=Cluster, y=Mean, fill=count_cluster)) +
#   geom_boxplot() +
#   #geom_jitter(color="black", size=0.4, alpha=0.9) + # si añado los puntos no se ve nada 
#   scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
#   theme(
#     legend.text = element_text(size=10),
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle("Average ranking of each GENE by cluster") +
#   xlab("Cluster") + 
#   scale_y_reverse() + theme_minimal() + labs(fill = "Gene Number")
# 


data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


##### plots 

#colors_palette <- c("#9770b1","#fdb6a1","#d387af","#666668","#fdfdd9")
#https://coolors.co/fec0aa-ec4e20-84732b-574f2a-1c3a13
colors_palette <- c("#FE5D26","#F2C078","#FAEDCA","#C1DBB3", "#7EBC89")
colors_palette <- c("#7EBC89","#C1DBB3","#FAEDCA","#F2C078", "#FE5D26")

p<-ggplot(mean_genes, aes(x=Cluster, y=Mean, fill=count_cluster)) +
  scale_fill_manual(values=colors_palette)+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") + 
  #ggtitle("Average ranking position of each gene") +
  geom_violin(trim=FALSE) + scale_y_reverse() + #theme_classic(base_size = 10) + 
  theme_bw(base_size = 15) + 
  labs(fill = "Number of genes") + 
  ylab("Gene Mean Ranking")  +
  stat_summary(fun.data=data_summary) + 
  theme(legend.position = c(0.79, 0.838), axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), plot.margin = margin(0.1,1,0.5,4, "cm"),
        legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))

p #+ stat_summary(fun.data=data_summary)


colors_palette <- c("#7EBC89","#C1DBB3","#FAEDCA","#F2C078", "#FE5D26")

gene_types$Var1 <- factor(gene_types$Var1 , levels = c( "pseudogenes", "ncRNA", "protein coding genes" ))
gene_types$cluster <- factor(gene_types$cluster)


# Graph
g <- ggplot(gene_types, aes(fill=cluster, y=Freq, x=Var1)) + 
  scale_fill_manual(values=colors_palette)+ theme_bw(base_size = 11) +
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~cluster, ncol=5) + coord_flip() + #theme( strip.background = element_blank() ) +
  theme(legend.position="none", axis.title.x = element_text(size = 20)) +   # theme_minimal() +
  xlab("") + ylab("") + #ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        plot.margin = margin(0,1,0,1, "cm"))
g


library("ggpubr")
ggarrange(p + font("x.text", size = 15) + font("y.text", size = 15), 
          g + font("x.text", size = 15) + font("y.text", size = 15), 
          heights = c(2, 0.4),
          ncol = 1, nrow = 2)

#p+geom_boxplot(width=0.1, fill="white") 
#p+geom_boxplot(width=0.08) 

##############3 voy a hacerlo también con ensamble (con BiomaRt): 

# Connect to the Ensembl database
library(biomaRt)
library(org.Hs.eg.db)  # For human genes; change depending on species
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_list <- rownames(disease_matrix)

# Retrieve information about gene type (coding, non-coding)
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list,
                         mart = ensembl)

length(gene_list)
nrow(gene_annotation)
head(gene_annotation)

# https://www.ensembl.org/info/genome/genebuild/biotypes.html

# Mapping table: Associate each biotype with its category: https://www.ensembl.org/info/genome/genebuild/biotypes.html
mapping <- data.frame(
  Var1 = c("IG_C_gene", "IG_C_pseudogene", "IG_J_gene", "IG_V_gene", "IG_V_pseudogene", 
           "lncRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "processed_pseudogene", 
           "protein_coding", "ribozyme", "rRNA", "rRNA_pseudogene", "scaRNA", "snoRNA", 
           "snRNA", "TEC", "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_V_pseudogene", 
           "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", 
           "transcribed_unprocessed_pseudogene", "unitary_pseudogene", 
           "unprocessed_pseudogene", "vault_RNA"),
  Category = c("IG C gene", "IG pseudogene", "IG J gene", "IG V gene", "IG pseudogene",
               "Long non-coding RNA (lncRNA)", "miRNA", "miscRNA", "MT-rRNA", "MT-tRNA", 
               "Processed pseudogene", "Protein coding", "Ribozyme", "rRNA", 
               "Pseudogene", "scaRNA", "snoRNA", "snRNA", "TEC", "TR C gene", 
               "TR J gene", "TR V gene", "TR pseudogene", "Transcribed pseudogene", 
               "Transcribed pseudogene", "Transcribed pseudogene", 
               "Unitary pseudogene", "Unprocessed pseudogene", "vaultRNA"),
  Label = c("Other", "Other", "Other", "Other", "Other", 
            "ncRNA", "ncRNA", "ncRNA", "Other", "Other", 
            "Pseudogene", "Protein coding", "ncRNA", "ncRNA", 
            "Pseudogene", "ncRNA", "ncRNA", "ncRNA", "Other", 
            "Other", "Other", "Other", "Other", 
            "Pseudogene", "Pseudogene", "Pseudogene", 
            "Pseudogene", "Pseudogene", "ncRNA")
)

dup_genes <- gene_annotation[duplicated(gene_annotation$hgnc_symbol),"hgnc_symbol"]
gene_annotation_dup_genes <- gene_annotation %>% dplyr::filter(hgnc_symbol %in% dup_genes)
length(unique(gene_annotation_dup_genes$hgnc_symbol))
#merge( gene_annotation_dup_genes , disease_matrix["Cluster"], by.x = "hgnc_symbol", by.y ="row.names" , all.x = TRUE, all.y = FALSE)

nrow(gene_annotation)
gene_annotation <- merge(gene_annotation, mapping,  by.x = "gene_biotype", by.y = "Var1", all.x = TRUE)
nrow(gene_annotation)

nrow(disease_matrix["Cluster"]) # 24757
gene_type_cluster <- merge(disease_matrix["Cluster"], gene_annotation, by.x = "row.names", by.y = "hgnc_symbol", all.x = TRUE, all.y = FALSE)
nrow(gene_type_cluster) # 24965

#

summary_gene_type_cluster <- gene_type_cluster %>%
  group_by(Cluster, Label) %>%
  summarise(Count = n(), .groups = "drop")

# los NAs son los que no se han anotado con Biomart: gene_annotation %>% dplyr::filter(hgnc_symbol == "AARS") ## este por ejemplo es un alias de AARS1

summary_gene_type_cluster$Cluster <- factor(summary_gene_type_cluster$Cluster)
summary_gene_type_cluster$Label <- factor(summary_gene_type_cluster$Label, levels = rev(c("Protein coding", "ncRNA", "Pseudogene", "Other")))

g <- summary_gene_type_cluster %>% dplyr::filter(Label == "Protein coding" | Label == "ncRNA" | Label == "Pseudogene")  %>% 
  ggplot(aes(fill=Cluster, y=Count, x=Label)) + 
  scale_fill_manual(values=colors_palette)+ theme_bw(base_size = 11) +
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~Cluster, ncol=5) + coord_flip() + #theme( strip.background = element_blank() ) +
  theme(legend.position="none", axis.title.x = element_text(size = 20)) +   # theme_minimal() +
  xlab("") + ylab("") + #ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        plot.margin = margin(0,1,0,1, "cm"))


library("ggpubr")
ggarrange(p + font("x.text", size = 15) + font("y.text", size = 15), 
          g + font("x.text", size = 15) + font("y.text", size = 15), 
          heights = c(2, 0.4),
          ncol = 1, nrow = 2)


### Algunos genes están anotados con varios gene_biotype distintos, supongo que porque tienen tránscritos alternativos. Vamos a ver, de los non-coding del cluster 1, cuales son únicamente non-coding
## y de los non-coding del cluster 5, cuales son también protein coding. 

genes_noncoding <- gene_type_cluster %>% dplyr::filter(Cluster == 1 & gene_biotype != "protein_coding")
gene_type_cluster %>% dplyr::filter(Row.names %in% genes_noncoding$Row.names)

genes_noncoding <- gene_type_cluster %>% dplyr::filter(Cluster == 5 & gene_biotype != "protein_coding")
gene_type_cluster %>% dplyr::filter(Row.names %in% genes_noncoding$Row.names) %>% dplyr::filter(gene_biotype == "protein_coding")

# a ver también que protein coding hay en el cluster 5

gene_type_cluster %>% dplyr::filter(Cluster == 5 & gene_biotype == "protein_coding")



###############3

stat_tops$cluster <- factor(stat_tops$cluster, levels = c_order$Cluster)

m <- ggplot(stat_tops, aes(x=median_no_0, y=Porc_panels_top)) + 
  #geom_point(aes(fill = factor(cluster)), color='#404040', shape=21, size=2)+ 
  geom_point(aes(colour = factor(cluster)))+ 
  geom_text_repel(data = filter(stat_tops, Porc_panels_top>0.75),aes(label=SYMBOL), size = 6) +
  scale_colour_manual(values=colors_palette)+
  #scale_fill_manual(values=c("#fdfdd9","#666668","#d387af","#fdb6a1","#9770b1"))+
  #scale_color_viridis(discrete = TRUE, alpha=0.9, option="A") +
  #ggtitle("Percentage of top position of each gene") + 
  labs(colour = "Cluster") + theme_light(base_size = 16) + guides(colour=guide_legend(ncol=5))+
  xlab("Gene Ranking Median") +  ylab("% Top Panel") + theme(legend.position = c(0.75, 0.85),
                                                             legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))

m 

m <- ggplot(stat_tops, aes(x=median_no_0, y=Porc_panels_top)) + 
  geom_point(aes(colour = factor(cluster)))+ 
  geom_text_repel(
    data = stat_tops %>% slice_max(Porc_panels_top, n = 10),
    aes(label = SYMBOL),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    force = 2
  ) +
  scale_colour_manual(values=colors_palette)+
  labs(colour = "Cluster") + 
  theme_light(base_size = 16) + 
  guides(colour=guide_legend(ncol=5))+
  xlab("Gene Ranking Median") +  
  ylab("% Top Panel") + 
  theme(
    legend.position = c(0.75, 0.85),
    legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979")
  )

m 

################
go_top <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_GO_Biological_Process_2025_table_TOP_Genes.txt", header = TRUE, sep = "\t", quote = "" )

go_top_005 <- go_top %>% dplyr::filter(Adjusted.P.value < 0.05)
go_top_005$NAME <- str_sub(go_top_005$Term,1,-14)
go_top_005[go_top_005$NAME == "Positive Regulation of Transcription by RNA Polymerase II" ,]$NAME <- "Positive Regulation of \nTranscription by RNA Polymerase II"
go_top_005[go_top_005$NAME == "Negative Regulation of Transcription by RNA Polymerase II" ,]$NAME <- "Negative Regulation of \nTranscription by RNA Polymerase II"
go_top_005[go_top_005$NAME == "Positive Regulation of Intracellular Signal Transduction" ,]$NAME <- "Positive Regulation of \nIntracellular Signal Transduction"


## vamos a redondear la leyenda
library(scales) # Required for scales::scientific

go_top_plot <- ggplot(go_top_005[c(1:15),], aes(x=reorder(NAME, Combined.Score), y=Combined.Score, color = Adjusted.P.value, size=Odds.Ratio)) +
  geom_point() +
  theme(panel.background = element_rect(fill='white', colour="black"),
        axis.text=element_text(size=14), axis.title=element_text(size=12), plot.title=element_text(size=12, face="bold", hjust=0.5),
        legend.key.size = unit(0.3, 'cm')) +
  scale_color_gradient(low="red", high="blue", 
                       labels = scales::scientific) +  # Correct usage of scales::scientific
  labs(x="", y="Combined.Score", size="Odds.Ratio", color="Adjusted.P.value") +
  coord_flip() + 
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0, "cm")
  )

go_top_plot


###################################3

stat_tops_toppgenes <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_stat_tops_vfinal.tsv", header = TRUE)
stat_tops_glowgenes <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", header = TRUE)

df_compare <- merge(
  stat_tops_toppgenes[, c("SYMBOL", "Porc_panels_top")],
  stat_tops_glowgenes[, c("SYMBOL", "Porc_panels_top")],
  by = "SYMBOL",
  suffixes = c("_topp", "_glow")
)

df_compare$diff <- df_compare$Porc_panels_top_topp - df_compare$Porc_panels_top_glow

head(df_compare[order(-abs(df_compare$diff)), ])


plot(df_compare$Porc_panels_top_topp,
     df_compare$Porc_panels_top_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)

abline(0, 1, col = "red")  # línea diagonal


###############################


#####

head(stat_tops)
sscore_glow <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", header = TRUE, sep = ",")
sscore_toppgenes <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_score_median.tsv", header = TRUE, sep = ",")

sscore_glow <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median_with_NAs.tsv", header = TRUE, sep = ",")
sscore_toppgenes <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_score_median_with_NAs.tsv", header = TRUE, sep = ",")


df_compare_sscore <- merge(
  sscore_toppgenes,
  sscore_glow,
  by = "Symbols",
  suffixes = c("_topp", "_glow")
)

df_compare_sscore$diff <- df_compare_sscore$Score_topp - df_compare_sscore$Score_glow
df_compare_sscore$diff_TOPP <- df_compare_sscore$median_ranking_topp - df_compare_sscore$min_ranking_topp
df_compare_sscore$diff_GLOW <- df_compare_sscore$median_ranking_glow - df_compare_sscore$min_ranking_glow


head(df_compare_sscore[order(-abs(df_compare_sscore$diff)), ])


plot(df_compare_sscore$Score_topp,
     df_compare_sscore$Score_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)

abline(0, 1, col = "red")  # línea diagonal

cor(df_compare_sscore$Score_topp, df_compare_sscore$Score_glow)

hist(df_compare_sscore$diff,
     breaks = 50,
     main = "Distribución de diferencias (ToppGene - Glowgenes)",
     xlab = "Diferencia en Porc_panels_top")
abline(v = 0, col = "red", lwd = 2)


#################3 vamos a ver qué son los genes esos extraños

genes_raros <- df_compare_sscore %>% dplyr::filter(diff_TOPP == 1655.0)
genes_raros <- df_compare_sscore %>% dplyr::filter(diff_GLOW == 821.0)


library(biomaRt)
library(org.Hs.eg.db)  # For human genes; change depending on species
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_list <- genes_raros$Symbols

# Retrieve information about gene type (coding, non-coding)
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list,
                         mart = ensembl)

length(gene_list)
nrow(gene_annotation)

as.data.frame(table(gene_annotation$gene_biotype))


############## vamos a filtrar solo los genes coding: 


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_list <- df_compare_sscore$Symbols

# Retrieve information about gene type (coding, non-coding)
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list,
                         mart = ensembl)



coding_genes <- gene_annotation %>% dplyr::filter(gene_biotype == "protein_coding")

df_compare_sscore_coding <- df_compare_sscore %>% dplyr::filter(Symbols %in% coding_genes$hgnc_symbol)

plot(df_compare_sscore_coding$Score_topp,
     df_compare_sscore_coding$Score_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)

cor(df_compare_sscore_coding$Score_topp, df_compare_sscore_coding$Score_glow)


plot(df_compare_sscore_coding$median_ranking_topp,
     df_compare_sscore_coding$median_ranking_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)

cor(df_compare_sscore_coding$median_ranking_topp, df_compare_sscore_coding$median_ranking_glow)


######## vamos a comparar las medianas

plot(df_compare_sscore$median_ranking_topp,
     df_compare_sscore$median_ranking_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)
abline(0, 1, col = "red")  # línea diagonal

cor(df_compare_sscore$median_ranking_topp, df_compare_sscore$median_ranking_glow)

## con el mínimo (el mejor valor) sale muy raro..

plot(df_compare_sscore$min_ranking_topp,
     df_compare_sscore$min_ranking_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)


median(
  as.numeric(disease_matrix["FAM181B", colnames(disease_matrix) != "Cluster"]),
  na.rm = TRUE
)

########################


df_compare_sscore$topp_diff_min_med <- df_compare_sscore$median_ranking_topp - df_compare_sscore$min_ranking_topp

df_compare_sscore_filtered <- df_compare_sscore%>% dplyr::filter(topp_diff_min_med != 1655)


plot(df_compare_sscore_filtered$Score_topp,
     df_compare_sscore_filtered$Score_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)

abline(0, 1, col = "red")  # línea diagonal

cor(df_compare_sscore_filtered$Score_topp, df_compare_sscore_filtered$Score_glow)



plot(df_compare_sscore_filtered$median_ranking_topp,
     df_compare_sscore_filtered$median_ranking_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)
abline(0, 1, col = "red")  # línea diagonal

cor(df_compare_sscore_filtered$median_ranking_topp, df_compare_sscore_filtered$median_ranking_glow)

## con el mínimo (el mejor valor) sale muy raro..

plot(df_compare_sscore_filtered$min_ranking_topp,
     df_compare_sscore_filtered$min_ranking_glow,
     xlab = "ToppGene",
     ylab = "Glowgenes",
     pch = 16)

###################
## vamos a ver los top y bottom genes según su score

sscore_glow[1:50,]
sscore_toppgenes[1:50,]


tail(sscore_glow, 50)
tail(sscore_toppgenes, 50)



# Obtener los 50 últimos genes de cada dataset
genes_glow <- tail(sscore_glow$Symbols, 50)
genes_topp <- tail(sscore_toppgenes$Symbols, 50)


intersect(genes_glow, genes_topp)


#install.packages("ggVennDiagram")
library(ggVennDiagram)

gene_lists <- list(
  GLOW = genes_glow,
  TOPPGENES = genes_topp
)

ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "blue", high = "red")



# Obtener los 50 primeros genes de cada dataset
genes_glow <- head(sscore_glow$Symbols, 50)
genes_topp <- head(sscore_toppgenes$Symbols, 50)


intersect(genes_glow, genes_topp)


#install.packages("ggVennDiagram")
library(ggVennDiagram)

gene_lists <- list(
  GLOW = genes_glow,
  TOPPGENES = genes_topp
)

ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "blue", high = "red")

### el primer percentil
# Calcular umbral del percentil 99
threshold_glow <- quantile(sscore_glow$Score, 0.99, na.rm = TRUE)

# Filtrar genes
genes_glow_top1 <- sscore_glow$Symbols[sscore_glow$Score >= threshold_glow]

threshold_topp <- quantile(sscore_toppgenes$Score, 0.99, na.rm = TRUE)

genes_topp_top1 <- sscore_toppgenes$Symbols[sscore_toppgenes$Score >= threshold_topp]


gene_lists <- list(
  GLOW = genes_glow_top1,
  TOPPGENES = genes_topp_top1
)

ggVennDiagram(gene_lists)

########
## el bottom 1% (percentil más bajo)
threshold_glow_low <- quantile(sscore_glow$Score, 0.01, na.rm = TRUE)

genes_glow_bottom1 <- sscore_glow$Symbols[sscore_glow$Score <= threshold_glow_low]

threshold_topp_low <- quantile(sscore_toppgenes$Score, 0.01, na.rm = TRUE)

genes_topp_bottom1 <- sscore_toppgenes$Symbols[sscore_toppgenes$Score <= threshold_topp_low]

length(genes_glow_bottom1)
length(genes_topp_bottom1)


gene_lists <- list(
  GLOW = genes_glow_bottom1,
  TOPPGENES = genes_topp_bottom1
)

ggVennDiagram(gene_lists)

################## CON LA MEDIANA ###########


# Obtener los 50 últimos genes de cada dataset por mediana
genes_topp <- tail(stat_tops[order(stat_tops$median, decreasing = TRUE),]$SYMBOL, 50)
genes_glow <- tail(stat_tops_GLOW[order(stat_tops_GLOW$median, decreasing = TRUE),]$SYMBOL, 50)


intersect(genes_glow, genes_topp)


#install.packages("ggVennDiagram")
library(ggVennDiagram)

gene_lists <- list(
  GLOW = genes_glow,
  TOPPGENES = genes_topp
)

ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "blue", high = "red")

# Obtener los 50 primeros genes de cada dataset por mediana
genes_topp <- head(stat_tops[order(stat_tops$median, decreasing = TRUE),]$SYMBOL, 50)
genes_glow <- head(stat_tops_GLOW[order(stat_tops_GLOW$median, decreasing = TRUE),]$SYMBOL, 50)


intersect(genes_glow, genes_topp)


#install.packages("ggVennDiagram")
library(ggVennDiagram)

gene_lists <- list(
  GLOW = genes_glow,
  TOPPGENES = genes_topp
)

ggVennDiagram(gene_lists) +
  scale_fill_gradient(low = "blue", high = "red")


###############3
# vamos a comparar los genes de los clusters

gene_topp_cluster1 <- stat_tops %>% dplyr::filter(cluster == "1") %>% dplyr::select(SYMBOL)
gene_glow_cluster1 <- stat_tops_GLOW %>% dplyr::filter(cluster == "1") %>% dplyr::select(SYMBOL)

gene_lists <- list(
  GLOW = gene_glow_cluster1$SYMBOL,
  TOPPGENES = gene_topp_cluster1$SYMBOL
)

ggVennDiagram(gene_lists)

###########
gene_topp_cluster5 <- stat_tops %>% dplyr::filter(cluster == "5") %>% dplyr::select(SYMBOL)
gene_glow_cluster5 <- stat_tops_GLOW %>% dplyr::filter(cluster == "5") %>% dplyr::select(SYMBOL)

gene_lists <- list(
  GLOW = gene_glow_cluster5$SYMBOL,
  TOPPGENES = gene_topp_cluster5$SYMBOL
)

ggVennDiagram(gene_lists)



