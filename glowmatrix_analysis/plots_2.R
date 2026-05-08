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

stat_tops  <- stat_tops %>%
  mutate(cluster = case_when(
    cluster == 4 ~ 1,
    cluster == 3 ~ 2,
    cluster == 0 ~ 3,
    cluster == 2 ~ 4,
    cluster == 1 ~ 5,
    TRUE ~ cluster # Default case (if needed)
  ))

gene_types  <- gene_types %>%
  mutate(cluster = case_when(
    cluster == 4 ~ 1,
    cluster == 3 ~ 2,
    cluster == 0 ~ 3,
    cluster == 2 ~ 4,
    cluster == 1 ~ 5,
    TRUE ~ cluster # Default case (if needed)
  ))


write.table(disease_matrix, "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_RENAMED_clustered_no_NA.tsv", row.names = FALSE, quote = FALSE, sep = "\t")


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

p + m
## Al final no hacemos este tsne, hacemos el de R. 
tsne["cluster"] <- stat_tops$cluster

t <- ggplot(tsne, aes(x=X0, y=X1)) + 
  geom_point(aes(colour = factor(cluster)))+ 
  scale_colour_manual(values=colors_palette)+
  #scale_fill_manual(values=c("#fdfdd9","#666668","#d387af","#fdb6a1","#9770b1"))+
  #scale_color_viridis(discrete = TRUE, alpha=0.9, option="A") +
  #ggtitle("t-SNE based on ranking position of each gene") + 
  labs(colour = "Cluster") +
  theme_minimal(base_size = 15) + 
  xlab("t-SNE1") +  ylab("t-SNE2") 
t

t + p / m + plot_layout(ncol = 2)
t / p + m + plot_layout(ncol = 2)

## Filtramos los genes



raw_glow_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv", header = TRUE)
rownames(raw_glow_matrix) <- raw_glow_matrix$SYMBOL


## vamos a quitar aquellos genes con NA cuyo mínimo valor sea mayor de 10.000. 

df_genes <- as.data.frame(raw_glow_matrix$SYMBOL)
colnames(df_genes) <- "SYMBOL"

#df_prueba <- as.data.frame(apply(raw_glow_matrix, 1, FUN = min,  na.rm = TRUE))

df_genes["min_value"] <- apply(raw_glow_matrix, 1, FUN = min,  na.rm = TRUE)
df_genes["NA_number"] <- rowSums(is.na(raw_glow_matrix))
df_genes["Cluster"] <- disease_matrix$Cluster

head(df_genes)
nrow(df_genes) # 24757

genes_entrezid <- mapIds(org.Hs.eg.db, rownames(disease_matrix), 'ENTREZID', 'SYMBOL')
genes_entrezid_df <- as.data.frame(genes_entrezid)
genes_entrezid_df["SYMBOL"] <- rownames(genes_entrezid_df)
df_genes <- merge(df_genes, genes_entrezid_df, by = "SYMBOL")

head(df_genes)
nrow(df_genes) # 24757

sum(is.na(df_genes$genes_entrezid))  # 884 NAs

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
attributes[1:5,]
attributes[attributes$page == "feature_page",]

gene_type <- getBM(attributes = c('entrezgene_id', 'gene_biotype'), 
                   filters = 'entrezgene_id', 
                   values = df_genes$genes_entrezid, 
                   mart = ensembl)


head(gene_type$entrezgene_id)
df_genes <- merge(df_genes, gene_type, by.x='genes_entrezid', by.y = "entrezgene_id", all.x = TRUE)
head(df_genes)
nrow(df_genes) # 24817
length(df_genes$SYMBOL)
length(unique(df_genes$SYMBOL)) # 24757

unique(df_genes$gene_biotype)
df_genes%>% dplyr::count(gene_biotype)
df_genes %>% dplyr::filter(is.na(gene_biotype))

filtered_genes <- df_genes %>% dplyr::filter(gene_biotype == "protein_coding" | min_value < 10000 | NA_number == "0")
length(unique(filtered_genes$SYMBOL)) # 20717  hemos quitado 4040 genes. 


nrow(filtered_genes)
filtered_genes %>% dplyr::count(Cluster)
unique(filtered_genes[c("SYMBOL", "Cluster")]) %>% dplyr::count(Cluster)
unique(df_genes[c("SYMBOL", "Cluster")]) %>% dplyr::count(Cluster)

filtered_genes <- merge(filtered_genes, stat_tops, by = "SYMBOL", all.x = TRUE)
filtered_genes <- filtered_genes[order(filtered_genes$Porc_panels_top, filtered_genes$median, decreasing = TRUE),] # ordenamos por top y mediana

#write.table(filtered_genes[filtered_genes["Cluster"] == 0, "SYMBOL"],
#            file = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/gene_sets/genes_clusterbottom.tsv", 
#            row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)


#write.table(filtered_genes[filtered_genes["Cluster"] == cluster_i, "SYMBOL"],
#            file = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/gene_sets/genes_clustertop.tsv", 
#            row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)

#

#go_bottom <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/enrichment/GO_Biological_Process_2023_table_bottom.txt", header = TRUE, sep = "\t", quote = "" )
go_bottom <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/enrichment/BOTTOM_allgenes_GO_Biological_Process_2025_table.txt", header = TRUE, sep = "\t", quote = "" )
go_top <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/enrichment/TOP_allgenes_GO_Biological_Process_2025_table.txt", header = TRUE, sep = "\t", quote = "" )

#reactome_bottom <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/enrichment/Reactome_2022_table_bottom.txt", header = TRUE, sep = "\t", quote = "" )
#reactome_top <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/enrichment/Reactome_2022_table_top.txt", header = TRUE, sep = "\t", quote = "" )

go_bottom_005 <- go_bottom %>% dplyr::filter(Adjusted.P.value < 0.05)

go_bottom_005$NAME <- str_sub(go_bottom_005$Term,1,-14)
go_bottom_005[go_bottom_005$NAME == "Detection of Chemical Stimulus Involved in Sensory Perception of Smell",]$NAME <- "Detection Of Chemical Stimulus \nInvolved In Sensory Perception Of Smell"
go_bottom_005[go_bottom_005$NAME == "Detection of Chemical Stimulus Involved in Sensory Perception" ,]$NAME <- "Detection Of Chemical Stimulus \nInvolved In Sensory Perception" 

go_bottom_plot <- ggplot(go_bottom_005, aes(x=reorder(NAME, Combined.Score), y=Combined.Score, color = Adjusted.P.value, size=Odds.Ratio)) +
  geom_point() +
  theme(panel.background = element_rect(fill='white', colour="black"),
        axis.text=element_text(size=14), axis.title=element_text(size=12), plot.title=element_text(size=12, face="bold", hjust=0.5),
        legend.key.size = unit(0.3, 'cm')) +
  scale_color_gradient(low="red", high="blue") +  # Correct usage of scales::scientific
  labs(x="", y="Combined.Score", size="Odds.Ratio", color="Adjusted.P.value") +
  coord_flip() + 
  theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0, "cm")
  )

go_bottom_plot

# reactome_bottom_005 <- reactome_bottom %>% dplyr::filter(Adjusted.P.value < 0.05)
# reactome_bottom_005$NAME <- sub(" R-HSA-\\d+$", "", reactome_bottom_005$Term)
# 
# reactome_bottom_005[reactome_bottom_005$NAME == "Expression And Translocation Of Olfactory Receptors",]$NAME <- "Expression And Translocation \n Of Olfactory Receptors"
# 
# react_bottom_plot <- ggplot(reactome_bottom_005, aes(x=reorder(NAME, Combined.Score), y=Combined.Score, color = Adjusted.P.value, size=Odds.Ratio)) +
#   geom_point() +
#   theme(panel.background = element_rect(fill='white', colour="black"), #legend.position="none", legend.position=c(-1,0.9),
#         axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
#         legend.key.size = unit(0.3, 'cm')) +
#   scale_color_gradient(low="red", high="blue") +  labs(x="",y="Combined.Score",  size="Odds.Ratio", color="Adjusted.P.value") +
#   coord_flip() + ggtitle("Reactome enrichment BOTTOM genes")  + theme(plot.margin = margin(1, 1, 1, 1, "cm"))



go_top_005 <- go_top %>% dplyr::filter(Adjusted.P.value < 0.05)
go_top_005$NAME <- str_sub(go_top_005$Term,1,-14)
go_top_005[go_top_005$NAME == "Positive Regulation of Transcription by RNA Polymerase II" ,]$NAME <- "Positive Regulation of \nTranscription by RNA Polymerase II"
go_top_005[go_top_005$NAME == "Negative Regulation of Transcription by RNA Polymerase II" ,]$NAME <- "Negative Regulation of \nTranscription by RNA Polymerase II"
go_top_005[go_top_005$NAME == "Positive Regulation of Intracellular Signal Transduction" ,]$NAME <- "Positive Regulation of \nIntracellular Signal Transduction"



# go_top_plot <- ggplot(go_top_005[c(1:15),], aes(x=reorder(NAME, Combined.Score), y=Combined.Score, color = Adjusted.P.value, size=Odds.Ratio)) +
#   geom_point() +
#   theme(panel.background = element_rect(fill='white', colour="black"), #legend.position="bottom" ,# legend.position=c(-1,0.9),
#         axis.text=element_text(size=12), axis.title=element_text(size=12), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
#         legend.key.size = unit(0.3, 'cm')) +
#   scale_color_gradient(low="red", high="blue") +  labs(x="",y="Combined.Score",  size="Odds.Ratio", color="Adjusted.P.value") +
#   coord_flip() + 
#   #ggtitle("GO-BP enrichment TOP genes") + 
#   theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) + 
#   theme(
#     legend.title = element_text(size = 12),       # Change the legend title font size
#     legend.text = element_text(size = 10),        # Change the legend text font size
#     legend.key.size = unit(0.5, "cm"),           # Change the size of the legend keys
#     legend.margin = margin(0, 0, 0, 0, "cm")     # Adjust the margin around the legend
#   )

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

# reactome_top_005 <- reactome_top %>% dplyr::filter(Adjusted.P.value < 0.05)
# reactome_top_005$NAME <- sub(" R-HSA-\\d+$", "", reactome_top_005$Term)
# 
# react_top_plot <- ggplot(reactome_top_005[c(1:20),], aes(x=reorder(NAME, Combined.Score), y=Combined.Score, color = Adjusted.P.value, size=Odds.Ratio)) +
#   geom_point() +
#   theme(panel.background = element_rect(fill='white', colour="black"), legend.position="bottom",# legend.position=c(-1,0.9),
#         axis.text=element_text(size=10), axis.title=element_text(size=10), plot.title=element_text(size=12, face="bold",hjust = 0.5 ),
#         legend.key.size = unit(0.3, 'cm')) +
#   scale_color_gradient(low="red", high="blue") +  labs(x="",y="Combined.Score",  size="Odds.Ratio", color="Adjusted.P.value") +
#   coord_flip() + ggtitle("Reactome enrichment TOP genes") + theme(plot.margin = margin(1, 1, 1, 1, "cm")) 
# 
# react_top_plot



pegos <- go_top_plot + go_bottom_plot
#preactomes <- react_top_plot + react_bottom_plot

grid.arrange(go_top_plot, go_bottom_plot, ncol = 2, heights = c(1), widths = c(1, 1))

pegos
#preactomes

t / p + m + pegos + plot_layout(ncol = 2) + plot_annotation(tag_levels = "I")

t / p + m + grid.arrange(go_top_plot , go_bottom_plot)  + plot_layout(ncol = 2) + plot_annotation(tag_levels = "I")


grid.arrange(go_top_plot , go_bottom_plot )

## cuadrar los tamaños del plot



ggdraw() +
  draw_plot(t, 0, .5, .5, .5) +
  draw_plot(p, 0, 0, .5, .5) +
  draw_plot(m, .5, .5, .5, .5) +
  draw_plot(pegos, .5, 0, .5, .5)
# draw_plot_label(c("A", "B", "C", "D"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)


rbind(c(1,1,2,2), c(1,1,2,2), c(3,3,4,5), c(3,3,4,5))
grid.arrange(t, p, m,  go_top_plot, go_bottom_plot, ncol = 2, 
             layout_matrix = rbind(c(1,1,2,2), c(1,1,2,2), c(3,3,4,5), c(3,3,4,5)))

rbind(c(1,1,1,2,2,2), c(1,1,1,2,2,2), c(3,3,4,4,5,5), c(3,3,4,4,5,5))
grid.arrange(t, m, p,  go_top_plot, go_bottom_plot, ncol = 2, 
             layout_matrix = rbind(c(1,1,1,2,2,2), c(1,1,1,2,2,2), c(3,3,4,4,5,5), c(3,3,4,4,5,5)))




#####

head(stat_tops)
sscore <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", header = TRUE, sep = ",")

sscore_porc <- merge(stat_tops, sscore, by.x = "SYMBOL", by.y = "Symbols")

m <- ggplot(sscore_porc, aes(x=Score, y=Porc_panels_top)) + 
  #geom_point(aes(fill = factor(cluster)), color='#404040', shape=21, size=2)+ 
  geom_point(aes(colour = factor(cluster)))+ 
  #geom_text_repel(data = filter(stat_tops, Porc_panels_top>0.75),aes(label=SYMBOL), size = 4) +
  scale_colour_manual(values=colors_palette)+
  #scale_fill_manual(values=c("#fdfdd9","#666668","#d387af","#fdb6a1","#9770b1"))+
  #scale_color_viridis(discrete = TRUE, alpha=0.9, option="A") +
  #ggtitle("Percentage of top position of each gene") + 
  labs(colour = "Cluster") + 
  theme_light(base_size = 16) + guides(colour=guide_legend(ncol=5))+
  #xlab("Gene Ranking Median") +  ylab("% Top Panel") + 
  theme(legend.position = c(0.75, 0.85),legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))

m 


###################################3

## exploración de los paneles de panelapp

glowmatrix <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = T, row.names = 1, stringsAsFactors = F, quote = "", check.names=F)
remove_cols <- c("COVID-19_research_GA", "Viral_resistance_GA", "Cluster")
glowmatrix <- glowmatrix[ , !(names(glowmatrix) %in% remove_cols)]

colnames(glowmatrix) <- gsub( "_GA$", "", colnames(glowmatrix))
head(glowmatrix)
nrow(glowmatrix)
ncol(glowmatrix)

panels_classification <- read.delim("~/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/panels/final_reclasified_panels.tsv", header = F, stringsAsFactors = F, quote = "", check.names=F)
colnames(panels_classification) <- c("classification","panel")

head(panels_classification)

# Find rows with 0 for each column
zero_rows <- apply(glowmatrix, 2, function(column) which(column == 0))
zero_counts <- apply(glowmatrix, 2, function(column) sum(column == 0))
zero_counts <- as.data.frame(zero_counts)

zero_counts <- merge(zero_counts, panels_classification, by.x = "row.names", by.y = "panel")
zero_counts <- zero_counts %>% group_by(classification) %>%
  mutate(median_zero_counts = median(zero_counts)) %>%
  ungroup()

ggplot(zero_counts, aes(x = reorder(classification, -median_zero_counts), y = log(zero_counts))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

View(zero_counts)

# number of panels by disease class
n_panels <- as.data.frame(table(panels_classification$classification))

ggplot(n_panels, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = NULL,
       y = NULL,
       fill = "Gene Type") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


zero_counts <- merge(zero_counts, n_panels, by.x = "classification", by.y = "Var1")
# Ensure the minimum and maximum zero_counts are included in the breaks
min_val <- min(zero_counts$zero_counts)
max_val <- max(zero_counts$zero_counts)

ggplot(zero_counts, aes(x = reorder(classification, -median_zero_counts), y = log2(zero_counts))) +
  geom_boxplot(aes(fill = Freq)) +
  scale_fill_gradient2(low = "#f7f0e7", high = "#389a3b", name = "Number of panels") +
  scale_y_continuous(
    name = "Number of genes",
    breaks = c(log2(c(min_val, max_val, 10, 50, 100)), scales::log_breaks(base = 2)(range(zero_counts$zero_counts))),
    labels = function(x) round(2^x, 1) # Convert log values back to original scale
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) +
  theme(
    legend.position = "bottom" # Change this to "bottom", "left", "right", or c(x, y)
  ) +
  labs(
    title = "Number of genes by panel",
    x = ""
  )

## ahora vamos a sacar cuántos paneles por gen hay: 

# Find rows with 0 for each row
zero_cols <- apply(glowmatrix, 1, function(row) which(row == 0))
zero_cols_counts <- apply(glowmatrix, 1, function(row) sum(row == 0))
zero_cols_counts <- as.data.frame(zero_cols_counts)
zero_cols_counts <- zero_cols_counts %>% dplyr::filter(zero_cols_counts>0)

df_zero_cols_counts <- as.data.frame(table(zero_cols_counts$zero_cols_counts))

ggplot(zero_cols_counts, aes(x = zero_cols_counts)) +
  geom_bar() +
  theme(
    legend.position = "bottom" # Change this to "bottom", "left", "right", or c(x, y)
  ) +
  labs(
    title = "Number of panel by gene",
    x = "Number of gene panels",
    y = "Number of genes"
  )




# Create a grouped variable: values > 20 become "20+"
zero_cols_counts <- zero_cols_counts %>%
  mutate(grouped = ifelse(zero_cols_counts >= 20, "20+", as.character(zero_cols_counts)))

# Explicit ordering: 1 → 20, then "20+"
zero_cols_counts$grouped <- factor(
  zero_cols_counts$grouped,
  levels = rev(c(as.character(1:19), "20+"))
)

# Plot
ggplot(zero_cols_counts, aes(x = grouped)) +
  geom_bar(fill = "#f0a773") +
  theme_minimal() +
  labs(
    title = "Number of gene panels by gene",
    x = "Number of gene panels",
    y = "Number of genes"
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))








