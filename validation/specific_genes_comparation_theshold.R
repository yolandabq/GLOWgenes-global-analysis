library(dplyr)

files <- Sys.glob("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/specific_genes/*.txt")

#files <- files[!grepl("bottom_0.75", files)]

genes_list <- list()

files


for (k in files){
  
  specific_genes_k <- read.table(k, sep = ",", header = TRUE)
  id <- basename(tools::file_path_sans_ext(k))
  id <- str_remove(id, "genes_panel_specific_")
  
  genes_list[[id]] <- specific_genes_k$SYMBOL
  
}

## intersección entre los 4 grupos
Reduce(intersect, genes_list)

### tamaño de los grupos 

sapply(genes_list, length)

###### matriz de solapamiento

overlap_matrix <- outer(names(genes_list), names(genes_list), Vectorize(function(a, b){
  length(intersect(genes_list[[a]], genes_list[[b]]))
}))

rownames(overlap_matrix) <- colnames(overlap_matrix) <- names(genes_list)

overlap_matrix


## proporción de solapamiento

prop_overlap <- outer(names(genes_list), names(genes_list), Vectorize(function(a, b){
  length(intersect(genes_list[[a]], genes_list[[b]])) / length(genes_list[[a]])
}))

rownames(prop_overlap) <- colnames(prop_overlap) <- names(genes_list)

prop_overlap

#### genes consistentes

library(dplyr)

gene_freq <- table(unlist(genes_list))

gene_freq <- as.data.frame(gene_freq)
colnames(gene_freq) <- c("gene", "n_lists")

table(gene_freq$n_lists)

##########3 qué proporción de genes aparece en más de 3 listas??


score <- sapply(names(genes_list), function(n){
  genes <- genes_list[[n]]
  mean(gene_freq$n_lists[match(genes, gene_freq$gene)] >= 3)
})

score

barplot(score, las = 2, main = "Proporción de genes compartidos")

#############333
## % de genes únicos y % de genes compartidos

genes_list <- list()

for (k in files){
  
  specific_genes_k <- read.table(k, sep = ",", header = TRUE)
  id <- basename(tools::file_path_sans_ext(k))
  id <- str_remove(id, "genes_panel_specific_")
  
  genes_list[[id]] <- specific_genes_k$SYMBOL
  
}


p <- as.data.frame(sapply(genes_list, length))
p 
#genes_list <- genes_list[!grepl("bottom_0.75", names(genes_list))]
#genes_list <- genes_list[!grepl("top_0.01_bottom_0.75", names(genes_list))]
#genes_list <- genes_list[!grepl("top_0.05_bottom_0.5", names(genes_list))]
genes_list <- genes_list[grepl("*_0.5", names(genes_list))]


gene_freq <- table(unlist(genes_list))


results <- sapply(names(genes_list), function(n){
  
  genes <- genes_list[[n]]
  freq <- gene_freq[genes]
  
  total <- length(genes)
  
  unique_pct <- sum(freq == 1) / total
  shared_pct <- sum(freq > 1) / total
  
  c(
    total_genes = total,
    pct_unique = unique_pct,
    pct_shared = shared_pct
  )
})

results <- t(results)
results


## genes compartidos en más de dos listas

results <- sapply(names(genes_list), function(n){
  
  genes <- genes_list[[n]]
  freq <- gene_freq[genes]
  
  total <- length(genes)
  
  c(
    total_genes = total,
    pct_unique = sum(freq == 1) / total,
    pct_in_2 = sum(freq == 2) / total,
    pct_in_3plus = sum(freq >= 3) / total
  )
})

t(results)



#########################33

library(dplyr)
library(tibble)

genes_list <- list()

for (k in files){
  
  specific_genes_k <- read.table(k, sep = ",", header = TRUE)
  id <- basename(tools::file_path_sans_ext(k))
  id <- str_remove(id, "genes_panel_specific_")
  
  genes_list[[id]] <- specific_genes_k$SYMBOL
  
}

genes_list <- genes_list[grepl("*_0.5", names(genes_list))]
genes_list

df_genes <- bind_rows(
  lapply(names(genes_list), function(g){
    tibble(
      gene = genes_list[[g]],
      group = g
    )
  })
)

diseases_textmining <- read.table("/home/yolanda/Downloads/human_disease_textmining_filtered.tsv", sep = "\t")

disease_tm_count <- diseases_textmining %>%
  group_by(V2) %>%                           # agrupa por gen
  summarise(n_diseases_TM = n_distinct(V3)) %>% # cuenta enfermedades únicas por DOID
  arrange(desc(n_diseases_TM))                  # opcional: ordenar de más a menos

df_genes <- merge(df_genes, disease_tm_count, 
                  by.x = "gene", by.y = "V2", 
                  all.x = TRUE)

df_genes$n_diseases_TM[is.na(df_genes$n_diseases_TM)] <- 0

ggplot(df_genes, aes(x = group, y = n_diseases_TM)) +
  geom_boxplot() +
  #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
  labs(x = "Threshold",
       y = "N de asociaciones en DISEASES TM",
       title = "Distribución del N de asociaciones según threshold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")


df_genes %>% dplyr::filter(n_diseases_TM != 0) %>% 
  ggplot(aes(x = group, y = n_diseases_TM)) +
    geom_boxplot() +
    #  geom_jitter(width = 0.2, alpha = 0.5) +  # opcional: puntos individuales
    labs(x = "Threshold",
         y = "N de asociaciones en DISEASES TM",
         title = "Distribución del N de asociaciones según threshold") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")


df_summary <- df_genes %>%
  group_by(group) %>%
  summarise(
    n_genes = n(),
    
    mean_diseases = mean(n_diseases_TM, na.rm = TRUE),
    median_diseases = median(n_diseases_TM, na.rm = TRUE),
    sd_diseases = sd(n_diseases_TM, na.rm = TRUE),
    
    # 👉 nuevos
    n_genes_nonzero = sum(n_diseases_TM > 0, na.rm = TRUE),
    mean_diseases_nonzero = mean(n_diseases_TM[n_diseases_TM > 0], na.rm = TRUE),
    
    min_diseases_no0 = ifelse(
      all(n_diseases_TM == 0),
      0,
      min(n_diseases_TM[n_diseases_TM > 0], na.rm = TRUE)
    ),
    
    max_diseases = max(n_diseases_TM, na.rm = TRUE),
    
    .groups = "drop"
  )

df_summary

library(tidyr)

# Pasar a formato largo
df_long <- df_summary %>%
  dplyr::select(group, mean_diseases, median_diseases) %>%
  tidyr::pivot_longer(
    cols = c(mean_diseases, median_diseases),
    names_to = "metric",
    values_to = "value"
  )

# Plot con líneas
ggplot(df_long, aes(x = group, y = value, color = metric, group = metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "Media y mediana de enfermedades por grupo",
    x = "Grupo",
    y = "Número de enfermedades",
    color = "Métrica"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(df_summary, aes(x = group, y = mean_diseases, group = 1)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  labs(
    title = "Average number of illnesses per group",
    x = "Group",
    y = "Mean diseases"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


############


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

mean(zero_counts$zero_counts)
median(zero_counts$zero_counts)

############

df <- imap_dfr(zero_rows, ~ tibble(
  panel = .y,
  gene = names(.x),
  value = as.numeric(.x)
))

df

df <- merge(df, panels_classification, by = "panel")


genes_per_class <- df %>%
                    group_by(classification) %>%
                    summarise(n_genes_unicos = n_distinct(gene)) %>%
                    arrange(desc(n_genes_unicos))


mean(genes_per_class$n_genes_unicos)
median(genes_per_class$n_genes_unicos)

n_genes <- length(unique(df$gene))


##########33

panels_list_tumour_cancer <- read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list.txt", header = F, stringsAsFactors = F, quote = "", check.names=F)
colnames(panels_list_tumour_cancer) <- c("Old_classification", "panel")
head(panels_list_tumour_cancer)
df <- merge(df, panels_list_tumour_cancer, by = "panel")


# Genes únicos por clasificación
cancer_genes <- df %>%
  filter(Old_classification == "Cancer Programme") %>%
  pull(gene) %>%
  unique()

tumour_genes <- df %>%
  filter(Old_classification == "Tumour syndromes") %>%
  pull(gene) %>%
  unique()

# Genes compartidos
shared_genes <- intersect(cancer_genes, tumour_genes)

shared_genes
length(shared_genes)

## genes totales
total_genes <- unique(c(cancer_genes, tumour_genes))
length(total_genes)

porc_shared <- length(shared_genes)/length(total_genes)

##################################

gene_list <- list(
  "Cancer Programme" = cancer_genes,
  "Tumour syndromes" = tumour_genes
)

# Venn
ggVennDiagram(gene_list) +
  scale_fill_gradient(low = "white", high = "steelblue")


##############3

genes_by_panel <- data.frame(table(df$panel))
