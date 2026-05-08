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


#stat_tops <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", header = TRUE)
disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA_7_clusters.tsv", header = TRUE)

rownames(disease_matrix) <- disease_matrix$SYMBOL
#tsne <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/tsne.tsv", header = TRUE)
#
#gene_types <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/gene_types.csv", header = TRUE, sep = ",")

### rename CLUSTERS FOR PLOT: 
disease_matrix <- disease_matrix %>%
  mutate(Cluster = case_when(
    Cluster == 3 ~ 0,
    Cluster == 6 ~ 1,
    Cluster == 5 ~ 2,
    Cluster == 4 ~ 3,
    Cluster == 2 ~ 4,
    Cluster == 0 ~ 5,
    Cluster == 1 ~ 6,
    TRUE ~ Cluster # Default case (if needed)
  ))


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
  # scale_fill_manual(values=colors_palette)+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") + 
  #ggtitle("Average ranking position of each gene") +
  geom_violin(trim=FALSE) + scale_y_reverse() + #theme_classic(base_size = 10) + 
  theme_bw(base_size = 15) + 
  labs(fill = "Number of genes") + 
  ylab("Gene Mean Ranking")  +
  stat_summary(fun.data=data_summary) + 
  theme(legend.position = c(0.79, 0.7), axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), plot.margin = margin(0.1,1,0.5,4, "cm"),
        legend.background = element_rect(size=0.5, linetype="dotdash", colour ="#797979"))

p #+ stat_summary(fun.data=data_summary)


############################################
## get gene type
########## con BiomaRt #####################
# Retrieve information about gene type (coding, non-coding)
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
  #scale_fill_manual(values=colors_palette)+ 
  theme_bw(base_size = 11) +
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~Cluster, ncol=7) + coord_flip() + #theme( strip.background = element_blank() ) +
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


#####################333


go_bottom <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/enrichment/GO_Biological_Process_2025_table_BOTTOM_CLUSTER7.txt", header = TRUE, sep = "\t", quote = "" )
go_top <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/enrichment/GO_Biological_Process_2025_table_TOP_CLUSTER7.txt", header = TRUE, sep = "\t", quote = "" )

#reactome_bottom <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/enrichment/Reactome_2022_table_bottom.txt", header = TRUE, sep = "\t", quote = "" )
#reactome_top <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/enrichment/Reactome_2022_table_top.txt", header = TRUE, sep = "\t", quote = "" )

go_bottom_005 <- go_bottom %>% dplyr::filter(Adjusted.P.value < 0.05)
## no hay ninguna función con un p ajustado < 0.05

go_top_005 <- go_top %>% dplyr::filter(Adjusted.P.value < 0.05)
go_top_005$NAME <- str_sub(go_top_005$Term,1,-14)
# go_top_005[go_top_005$NAME == "Positive Regulation of Transcription by RNA Polymerase II" ,]$NAME <- "Positive Regulation of \nTranscription by RNA Polymerase II"
# go_top_005[go_top_005$NAME == "Negative Regulation of Transcription by RNA Polymerase II" ,]$NAME <- "Negative Regulation of \nTranscription by RNA Polymerase II"
# go_top_005[go_top_005$NAME == "Positive Regulation of Intracellular Signal Transduction" ,]$NAME <- "Positive Regulation of \nIntracellular Signal Transduction"

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

















