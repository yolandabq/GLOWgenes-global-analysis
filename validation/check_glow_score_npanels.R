GAR_genes <- list.files(
  path = "/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/panels",
  pattern = "_GAR.csv",
  recursive = TRUE,
  full.names = TRUE
)

disease_classes <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/GLOWgenes_RDs/panels/final_reclasified_panels.tsv", sep = "\t")

# Obtener solo nombres de archivo
files <- basename(GAR_genes)
# Quitar el sufijo _GAR.csv
names_clean <- sub("_GAR\\.csv$", "", files)
names_clean

GAR_genes <- as.data.frame(GAR_genes)
GAR_genes$panel <- names_clean

GAR_genes <- GAR_genes %>% dplyr::filter(panel %in% disease_classes$V2)


for (g in GAR_genes$GAR_genes){
  
  print(g)
  df_panelApp <- read.table(
    file = g,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  red_genes <- df_panelApp[grepl("Red", df_panelApp$evidence), ]
  
  print(nrow(red_genes))
  
}





df_red_all <- data.frame(
  gene = character(),
  panel = character(),
  stringsAsFactors = FALSE
)

for (g in GAR_genes$GAR_genes){
  
  print(g)
  
  df_panelApp <- read.table(
    file = g,
    header = TRUE,
    sep = "\t",
    fill = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  red_genes <- df_panelApp[grepl("Red", df_panelApp$evidence), ]
  
  # nombre del panel (archivo sin ruta ni sufijo)
  panel_name <- sub("_GAR\\.csv$", "", basename(g))
  
  if(nrow(red_genes) > 0){
    temp <- data.frame(
      gene = red_genes$entity_name,
      panel = panel_name,
      stringsAsFactors = FALSE
    )
    
    df_red_all <- rbind(df_red_all, temp)
  }
}

df_red_all

df_red_all <- merge(df_red_all, disease_classes, by.x = "panel", by.y = "V2")


gene_counts <- df_red_all %>%
  group_by(gene) %>%
  summarise(n_panels = n_distinct(panel)) %>%
  arrange(desc(n_panels))

gene_counts

## en vez de por panel, por clase de enfermedad (habrá paneles de la misma clase que pueden compartir genes red)
gene_counts <- df_red_all %>%
  group_by(gene) %>%
  summarise(n_panels = n_distinct(V1)) %>%
  arrange(desc(n_panels))

gene_counts



df_score = read.table("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", 
                      head=TRUE, sep = ",")


df_genes_score <- merge(gene_counts, df_score, by.x = "gene", by.y = "Symbols", all.x = TRUE)
df_genes_score$n_panels <- factor(df_genes_score$n_panels, levels = sort(unique(df_genes_score$n_panels)))*

## hay 82 genes que no tienen score (82 genes que no hemos pillado)
df_genes_score[!is.finite(df_genes_score$n_panels) | !is.finite(df_genes_score$Score), ]


ggplot(df_genes_score, aes(x=n_panels, y=Score, fill=n_panels)) +
  geom_boxplot(alpha = 0.5) + #geom_jitter() +
  theme_classic(base_size = 15) + labs(fill = "N panel") + guides(fill = guide_legend(ncol = 2))+ 
  labs(title="Gene Disease Specificity Score",x="Number of disease-classes panel", y = "Gene Disease Specificity Score")


df_genes_score %>% dplyr::filter(Score > 0.5) %>% 
  ggplot(aes(x=n_panels, y=Score, fill=n_panels)) +
  geom_boxplot(alpha = 0.5) + #geom_jitter() +
  theme_classic(base_size = 15) + labs(fill = "N panel") + guides(fill = guide_legend(ncol = 2))+ 
  labs(title="Gene Disease Specificity Score",x="Number of disease-classes panel", y = "Gene Disease Specificity Score")



as.data.frame(table(df_genes_score$n_panels))








