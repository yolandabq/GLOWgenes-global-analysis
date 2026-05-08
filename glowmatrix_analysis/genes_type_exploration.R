# Install packages if you haven't already
# if (!requireNamespace("biomaRt", quietly = TRUE)) {
#   install.packages("biomaRt")
# }
# if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#   BiocManager::install("org.Hs.eg.db")
# }
# if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
#   BiocManager::install("GenomicFeatures")
# }

# Load libraries
library(biomaRt)
library(org.Hs.eg.db)  # For human genes; change depending on species
library(GenomicFeatures)
library(ggplot2)
library(dplyr)
library(stringr)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Your gene list (replace with your actual gene names)
disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = TRUE)
gene_list <- disease_matrix$SYMBOL


########## con BiomaRt #####################

# Retrieve information about gene type (coding, non-coding)
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list,
                         mart = ensembl)

# View the results
print(gene_annotation)
as.data.frame(table(gene_annotation$gene_biotype))
sum(as.data.frame(table(gene_annotation$gene_biotype))$Freq) # aquí hay más información (24024), pero habría que mirarla. Quizás con la otra para el charplot ya esta bien. 
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


# Merge the original dataframe with the mapping to add the category
merged_df <- merge(as.data.frame(table(gene_annotation$gene_biotype)), mapping, by = "Var1", all.x = TRUE)
head(merged_df)

# Summarize frequencies by category
result <- aggregate(Freq ~ Label, data = merged_df, sum)

head(result)

ggplot(result, aes(x = "", y = Freq, fill = Label)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = NULL,
       y = NULL,
       fill = "Gene Type") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())


########## con AnnotationDbi #####################   
## sale un poco distinto. Con BiomaRt se anota más genes. 

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Get gene type annotation using Entrez IDs
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")

# View results
print(gene_info)

gene_type <- as.data.frame(table(gene_info$GENETYPE))
sum(gene_type$Freq) # 23859 genes que se han podido anotar con el tipo de gen
length(gene_list) # 24757 (los mismos que ne glowmatrix)

gene_list[! gene_list %in% gene_info$SYMBOL] # estos son los genes que no se han podido anotar

gene_info %>% dplyr::filter(SYMBOL == "C2orf78" )


# Plot the data
ggplot(gene_type, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = "Gene Type",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()


# Combine specific rows into "ncRNA"
combined_ncRNA <- gene_type %>%
  filter(Var1 %in% c("ncRNA", "rRNA", "scRNA", "snoRNA", "snRNA", "tRNA")) %>%
  summarise(
    Var1 = "ncRNA",
    Freq = sum(Freq)
  )

# Filter out the original rows that were grouped
gene_type <- gene_type %>%
  filter(!Var1 %in% c("ncRNA", "rRNA", "scRNA", "snoRNA", "snRNA", "tRNA"))

# Add the new combined "ncRNA" row back
gene_type <- bind_rows(gene_type, combined_ncRNA)

gene_type_plot <- gene_type %>% dplyr::filter(Var1 == "protein-coding" | Var1 == "pseudo" | Var1 == "ncRNA")
gene_type_plot$Var1 <- str_replace(gene_type_plot$Var1, "protein-coding", 'protein coding genes')
gene_type_plot$Var1 <- str_replace(gene_type_plot$Var1, "pseudo", 'pseudogenes')
gene_type_plot$Var1 <- str_replace(gene_type_plot$Var1, "ncRNA_total", 'ncRNA')

#gene_type_plot$Var1 <- factor(gene_type_plot$Var1 , levels = c("protein-coding", "ncRNA_total", "pseudo"))
gene_type_plot$Var1 <- factor(gene_type_plot$Var1 , levels = c( "pseudogenes", "ncRNA", "protein coding genes" ))


ggplot(gene_type_plot, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = NULL,
       y = NULL,
       fill = "Gene Type") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())



###### Genes de PanelApp


# Your gene list (replace with your actual gene names)
genes_panelapp <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/gene_sets/all_genes_panelapp.tsv", header = TRUE)
gene_list <- genes_panelapp$genes
length(gene_list) # 4414 unique genes in panelapp

########## con BiomaRt #####################   
# Retrieve information about gene type (coding, non-coding)
gene_annotation <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list,
                         mart = ensembl)

# View the results
print(gene_annotation)
as.data.frame(table(gene_annotation$gene_biotype))
sum(as.data.frame(table(gene_annotation$gene_biotype))$Freq)  # 4315 genes annotated

# Merge the original dataframe with the mapping to add the category
merged_df <- merge(as.data.frame(table(gene_annotation$gene_biotype)), mapping, by = "Var1", all.x = TRUE)
merged_df
# Summarize frequencies by category
result <- aggregate(Freq ~ Label, data = merged_df, sum)
result

ggplot(result, aes(x = "", y = Freq, fill = Label)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = NULL,
       y = NULL,
       fill = "Gene Type") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

########## con AnnotationDbi #####################   
# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Get gene type annotation using Entrez IDs
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")

# View results
print(gene_info)

gene_type <- as.data.frame(table(gene_info$GENETYPE))
sum(gene_type$Freq)
length(gene_list)
gene_type

# Plot the data
ggplot(gene_type, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Gene Type Frequency",
       x = "Gene Type",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()




##### genes by cluster



disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = TRUE)
gene_list0 <- disease_matrix %>% dplyr::filter(Cluster == 0) %>% dplyr::select(SYMBOL)
gene_list1 <- disease_matrix %>% dplyr::filter(Cluster == 1) %>% dplyr::select(SYMBOL)
gene_list2 <- disease_matrix %>% dplyr::filter(Cluster == 2) %>% dplyr::select(SYMBOL)
gene_list3 <- disease_matrix %>% dplyr::filter(Cluster == 3) %>% dplyr::select(SYMBOL)
gene_list4 <- disease_matrix %>% dplyr::filter(Cluster == 4) %>% dplyr::select(SYMBOL)

########## con BiomaRt #####################   
# Retrieve information about gene type (coding, non-coding)
gene_annotation0 <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list0,
                         mart = ensembl)

# Retrieve information about gene type (coding, non-coding)
gene_annotation1 <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list1,
                         mart = ensembl)

# Retrieve information about gene type (coding, non-coding)
gene_annotation2 <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list2,
                         mart = ensembl)

# Retrieve information about gene type (coding, non-coding)
gene_annotation3 <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list3,
                         mart = ensembl)

# Retrieve information about gene type (coding, non-coding)
gene_annotation4 <- getBM(attributes = c('hgnc_symbol', 'gene_biotype'),
                         filters = 'hgnc_symbol',
                         values = gene_list4,
                         mart = ensembl)
# View the results
as.data.frame(table(gene_annotation0$gene_biotype))
as.data.frame(table(gene_annotation1$gene_biotype))
as.data.frame(table(gene_annotation2$gene_biotype))
as.data.frame(table(gene_annotation3$gene_biotype))
as.data.frame(table(gene_annotation4$gene_biotype))

########## con AnnotationDbi #####################   
# Convert gene symbols to Entrez IDs
entrez_ids0 <- mapIds(org.Hs.eg.db, keys = gene_list0$SYMBOL, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids1 <- mapIds(org.Hs.eg.db, keys = gene_list1$SYMBOL, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids2 <- mapIds(org.Hs.eg.db, keys = gene_list2$SYMBOL, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids3 <- mapIds(org.Hs.eg.db, keys = gene_list3$SYMBOL, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids4 <- mapIds(org.Hs.eg.db, keys = gene_list4$SYMBOL, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")


# Get gene type annotation using Entrez IDs
gene_info0 <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids0, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")
gene_info1 <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids1, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")
gene_info2 <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids2, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")
gene_info3 <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids3, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")
gene_info4 <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids4, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID")




# View results

gene_type0 <- as.data.frame(table(gene_info0$GENETYPE))
sum(gene_type0$Freq)
length(gene_list0)
gene_type0
gene_type0$cluster <- 0


gene_type1 <- as.data.frame(table(gene_info1$GENETYPE))
sum(gene_type1$Freq)
length(gene_list1)
gene_type1



gene_type1$cluster <- 1


gene_type2 <- as.data.frame(table(gene_info2$GENETYPE))
sum(gene_type2$Freq)
length(gene_list2)
gene_type2
gene_type2$cluster <- 2


gene_type3 <- as.data.frame(table(gene_info3$GENETYPE))
sum(gene_type3$Freq)
length(gene_list3)
gene_type3
gene_type3$cluster <- 3


gene_type4 <- as.data.frame(table(gene_info4$GENETYPE))
sum(gene_type4$Freq)
length(gene_list4)
gene_type4
gene_type4$cluster <- 4


gene_type_all <- rbind(gene_type0, gene_type1, gene_type2, gene_type3, gene_type4)

nc_rna <- gene_type_all %>% dplyr::filter(Var1 == "ncRNA" | Var1 == "rRNA" | Var1 == "scRNA"
                                       | Var1 == "snoRNA" | Var1 == "snRNA" | Var1 == "tRNA") %>% 
  dplyr::group_by(cluster) %>%
  dplyr::summarise(sum(Freq))

nc_rna_df <- as.data.frame(nc_rna)
colnames(nc_rna_df) <- c("cluster", "Freq")
nc_rna_df$Var1 <- "ncRNA_total"

gene_type_all <- rbind(gene_type_all, nc_rna_df[c("Var1", "Freq", "cluster")])

gene_type_plot <- gene_type_all %>% dplyr::filter(Var1 == "protein-coding" | Var1 == "pseudo" | Var1 == "ncRNA_total")
gene_type_plot$Var1 <- str_replace(gene_type_plot$Var1, "protein-coding", 'protein coding genes')
gene_type_plot$Var1 <- str_replace(gene_type_plot$Var1, "pseudo", 'pseudogenes')
gene_type_plot$Var1 <- str_replace(gene_type_plot$Var1, "ncRNA_total", 'ncRNA')

#gene_type_plot$Var1 <- factor(gene_type_plot$Var1 , levels = c("protein-coding", "ncRNA_total", "pseudo"))
gene_type_plot$Var1 <- factor(gene_type_plot$Var1 , levels = c( "pseudogenes", "ncRNA", "protein coding genes" ))
gene_type_plot$cluster <- factor(gene_type_plot$cluster, levels = c(4,3,0,2,1) )

write.csv(gene_type_plot, "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/gene_types.csv", row.names=FALSE)

#gene_type_plot <- gene_type_plot %>% group_by(cluster) %>% mutate(percentage = (Freq / sum(Freq)) * 100) %>% ungroup()



colors_palette <- c("#7EBC89","#C1DBB3","#FAEDCA","#F2C078", "#FE5D26")


# Graph
g <- ggplot(gene_type_plot, aes(fill=cluster, y=Freq, x=Var1)) + 
  scale_fill_manual(values=colors_palette)+ theme_bw() +
  geom_bar(position="dodge", stat="identity")+
  facet_wrap(~cluster, ncol=5) + coord_flip() + #theme( strip.background = element_blank() ) +
  theme(legend.position="none") +   # theme_minimal() +
  xlab("") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
g



library("ggpubr")
ggarrange(p + font("x.text", size = 10) + font("y.text", size = 10), 
          g + font("x.text", size = 15) + font("y.text", size = 12), 
          heights = c(2, 0.5),
          ncol = 1, nrow = 2)

###


glow_panelappgenes<- disease_matrix %>% dplyr::filter(SYMBOL %in% genes_panelapp$genes) 
as.data.frame(table(glow_panelappgenes$Cluster))










