library(dplyr)

disease_matrix <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/glow_matrix_clustered_no_NA.tsv", header = TRUE)
rownames(disease_matrix) <- disease_matrix$SYMBOL
#

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


disease_matrix_c5 <- disease_matrix %>% dplyr::filter(Cluster == "5")
c5_genes <- disease_matrix_c5$SYMBOL

write.table(c5_genes,
            file = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/enrichment/genes_cluster_5.tsv", 
            row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)

disease_matrix_c1 <- disease_matrix %>% dplyr::filter(Cluster == "1")
c1_genes <- disease_matrix_c1$SYMBOL

write.table(c1_genes,
            file = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/enrichment/genes_cluster_1.tsv", 
            row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)
