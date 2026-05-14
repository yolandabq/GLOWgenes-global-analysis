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
rownames(disease_matrix) <- disease_matrix$SYMBOL
sscore_glow <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", header = TRUE, sep = ",")
stats <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", header = TRUE, sep = "\t")

#
#gene_types <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/gene_types.csv", header = TRUE, sep = ",")

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


######################3

# numero de gos asociados al cluster 1

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# numero de gos en los top
go_data <- getBM(
  attributes = c("hgnc_symbol", "go_id", "namespace_1003"),
  filters    = "hgnc_symbol",
  values     = disease_matrix$SYMBOL,
  mart       = mart
)

# Filtrar solo Biological Process
go_bp <- go_data %>%
  dplyr::filter(namespace_1003 == "biological_process", go_id != "") %>%
  dplyr::distinct()

##Contar nº de GO-BP por gen
library(dplyr)

go_bp_counts <- go_bp %>%
  group_by(hgnc_symbol) %>%
  summarise(n_GO_BP = n(), .groups = "drop")

go_bp_counts

mean(go_bp_counts$n_GO_BP)

go_bp_counts <- merge(go_bp_counts,sscore_glow , by.x = "hgnc_symbol", "Symbols")

genes_highlight <- c("AKT1", "GAPDH", "HRAS")

# Plot base
plot(go_bp_counts$n_GO_BP,
     go_bp_counts$Score,
     xlab = "Nº GOs",
     ylab = "SGDS",
     pch = 16)

cor(go_bp_counts$n_GO_BP, go_bp_counts$Score)


ggplot(go_bp_counts, aes(x = n_GO_BP , y = Score)) +
  geom_point(alpha = 0.3, colour = "red") +
  
  # Resaltar genes
  geom_point(
    data = go_bp_counts %>% filter(hgnc_symbol %in% genes_highlight),
    size = 2,
    colour = "green",
  ) +
  
  
  # Etiquetas
  geom_text_repel(
    data = go_bp_counts %>% filter(hgnc_symbol %in% genes_highlight),
    aes(label = hgnc_symbol),
    size = 4,
    show.legend = FALSE,
    colour = "black"
  ) 

go_bp_counts <- merge(go_bp_counts,stats , by.x = "hgnc_symbol", "SYMBOL")


ggplot(go_bp_counts, aes(x = n_GO_BP , y = Porc_panels_top)) +
  geom_point(alpha = 0.3, colour = "red") +
  
  # Resaltar genes
  geom_point(
    data = go_bp_counts %>% filter(hgnc_symbol %in% genes_highlight),
    size = 2,
    colour = "green",
  ) +
  
  
  # Etiquetas
  geom_text_repel(
    data = go_bp_counts %>% filter(hgnc_symbol %in% genes_highlight),
    aes(label = hgnc_symbol),
    size = 4,
    show.legend = FALSE,
    colour = "black"
  ) 




# Calcular los cuartiles
quantiles <- quantile(go_bp_counts$n_GO_BP, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
quantiles

# Crear variable categórica usando cut con los cuartiles
go_bp_counts$N_GOs_quartiles <- cut(go_bp_counts$n_GO_BP,
                                      breaks = quantiles,
                                      include.lowest = TRUE,
                                      labels = c("Q1","Q2","Q3","Q4"))

as.data.frame(table(go_bp_counts$N_GOs_quartiles))

ggplot(go_bp_counts, aes(x = N_GOs_quartiles, y = Score, fill = N_GOs_quartiles)) +
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

go_bp_counts$N_GOs_quartiles <- cut(
  go_bp_counts$n_GO_BP,
  breaks = quantiles,
  include.lowest = TRUE,
  labels = labels_quartiles
)

ggplot(go_bp_counts, aes(x = N_GOs_quartiles, y = Score, fill = N_GOs_quartiles)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, alpha = 0.1) +  # muestra puntos individuales
  labs(x = "N GOs",
       y = "SGDS",
       title = "Distribución del Score según número de GOs",
       fill = "Quantiles") +
  theme_minimal(base_size = 18) +
  scale_fill_brewer(palette = "Set3")

######################

sscore_glow <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", header = TRUE, sep = ",")
#betweenness <- read.table("~/tblab/yolanda/GLOWgenes/panelAPP/revision/bewteenness_specific_general.csv", header = TRUE, sep = ",")
betweenness <- read.table("~/tblab/yolanda/GLOWgenes/panelAPP/revision/bewteenness_specific_general_3_k50.csv", header = TRUE, sep = ",")

df_long <- betweenness %>%
  pivot_longer(
    cols = -c(Gene, type),
    names_to = "metric",
    values_to = "value"
  )

df_long[is.na(df_long)] <- 0

ggplot(df_long, aes(x = type, y = value, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  labs(
    title = "Comparación General vs Specific por métrica",
    x = "Tipo de gen",
    y = "Valor"
  )


df_long <- df_long %>%
  mutate(value_log = log10(value + 1))

library(stringr)

df_long <- df_long %>%
  mutate(
    metric = str_replace(
      metric,
      "^(coexpression|coregulation|drug|pathway|database|function|phenotype|coessent|complex|coocurrence|fitness|localization|genetic|fusion|physical|neighborhood|tfregulon|experimental|regulon|textmining)",
      "\\1-"
    )
  )

df_long <- df_long %>%
  mutate(
    metric = str_replace(
      metric,
      "AVANAwEXT",
      "AVANA"
    )
  )

df_long <- df_long %>%
  mutate(
    metric = str_replace(
      metric,
      "AVANAAEXT",
      "AVANA"
    )
  )


df_long <- df_long %>%
  mutate(
    metric = str_replace(
      metric,
      "CRISPRw",
      "CRISPR"
    )
  )


df_long <- df_long %>%
  mutate(
    type = str_replace(
      type,
      "general",
      "global-genes"
    )
  )

df_long <- df_long %>%
  mutate(
    type = str_replace(
      type,
      "specific",
      "specific-genes"
    )
  )





ggplot(df_long, aes(x = type, y = value_log, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() + theme(
    legend.position = "none"
  ) +
  labs(x = NULL)

df_long %>% 
  filter(metric %in% c(
    "coessent-AVANA",
    "fitness-AVANA",
    "function-GOBP",
    "textmining-STRING"
  )) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "function-GOBP",
        "textmining-STRING",
        "coessent-AVANA",
        "fitness-AVANA"
      )
    )
  ) %>%
  ggplot(aes(x = type, y = value_log, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free", nrow = 2) +
  labs(
    title = NULL,
    y = "Betweenness (log)",
    x = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none"
  )

head(sscore_glow)
head(df_long)

df_long <- merge(df_long, sscore_glow, by.x = "Gene", by.y = "Symbols", all.x=TRUE, all.y=FALSE)

ggplot(df_long, aes(x = Score, y = value, colour = type)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  labs(
    title = "Relación entre Score y valor por métrica",
    x = "SGDS",
    y = "Valor"
  )



###################################3
betweenness <- read.table("~/tblab/yolanda/GLOWgenes/panelAPP/revision/bewteenness_specific_general_3_k50.csv", header = TRUE, sep = ",")
betweenness <- read.table("~/tblab/yolanda/GLOWgenes/panelAPP/revision/bewteenness_HPO_ext_no_weight.csv", header = TRUE, sep = ",")

df_long


df_long <- betweenness %>%
  pivot_longer(
    cols = -c(Gene, type),
    names_to = "metric",
    values_to = "value"
  )

df_long[is.na(df_long)] <- 0

df_long <- merge(df_long, sscore_glow, by.x = "Gene", by.y = "Symbols", all.x=TRUE, all.y=FALSE)
df_long <- merge(df_long, stats, by.x = "Gene", by.y = "SYMBOL", all.x=TRUE, all.y=FALSE)



genes_highlight <- c("AKT1", "GAPDH", "HRAS")

g <- ggplot(df_long, aes(x = Score , y = value)) +
  geom_point(alpha = 0.3, colour = "red") +
  
  # Resaltar genes
  geom_point(
    data = df_long %>% filter(Gene %in% genes_highlight),
    size = 2,
    colour = "green",
  ) +
  
  # Resaltar genes
  geom_point(
    data = df_long %>% filter(type=="specific"),
    size = 2,
    colour = "blue", alpha = 0.4
  ) +
  
  # Etiquetas
  geom_text_repel(
    data = df_long %>% filter(Gene %in% genes_highlight),
    aes(label = Gene),
    size = 4,
    show.legend = FALSE,
    colour = "black"
  ) +
  
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  labs(
    title = "Relación entre Score y valor por métrica"
  )


ggsave("~/tblab/yolanda/GLOWgenes/panelAPP/revision/betweenness_plot.png", plot = g)

df_long %>%
  filter(metric == "coexpressionCOXPRESdbEXT") %>%
  arrange(desc(value))


# Filtrar métrica
df_metric <- df_long %>%
  filter(metric == "coexpressionCOXPRESdbEXT")

df_metric <- df_long %>%
  filter(metric == "coexpressionSTRING")

df_metric <- df_long %>%
  filter(metric == "fitnessAVANAAEXT")

df_metric <- df_long %>%
  filter(metric == "physicalHIunion")

# Percentil 90
p90 <- quantile(df_metric$value, 0.90, na.rm = TRUE)

# Genes en percentil 90
df_p90 <- df_metric %>%
  filter(value >= p90)

## genes specific en el percentil 90

df_p90 %>% dplyr::filter(type == "specific") %>% arrange(desc(value))
df_p90 %>% dplyr::filter(type == "general") %>% arrange(desc(value))

### qué porcentaje del percentil 90 es specific y general

df_p90 %>%
  count(type) %>%
  mutate(percentage = 100 * n / sum(n))


#### de todos los specific, cuantos estan en el percentil 90

n_specific_total <- df_metric %>%
  filter(type == "specific") %>%
  nrow()

n_specific_p90 <- df_p90 %>%
  filter(type == "specific") %>%
  nrow()

percentage_specific_in_p90 <- 100 * n_specific_p90 / n_specific_total

percentage_specific_in_p90

########### PROBAMOS CON EL DEGREE

degree <- read.table("~/tblab/yolanda/GLOWgenes/panelAPP/revision/degree_no_weight.csv", header = TRUE, sep = ",")
sscore_glow <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", header = TRUE, sep = ",")
stats <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", header = TRUE, sep = "\t")


df_long <- degree %>%
  pivot_longer(
    cols = -c(Gene, type),
    names_to = "metric",
    values_to = "value"
  )

df_long[is.na(df_long)] <- 0

df_long <- merge(df_long, sscore_glow, by.x = "Gene", by.y = "Symbols", all.x=TRUE, all.y=FALSE)
df_long <- merge(df_long, stats, by.x = "Gene", by.y = "SYMBOL", all.x=TRUE, all.y=FALSE)



genes_highlight <- c("AKT1", "GAPDH", "HRAS")

g <- ggplot(df_long, aes(x = Score , y = value)) +
  geom_point(alpha = 0.3, colour = "red") +
  
  # Resaltar genes
  geom_point(
    data = df_long %>% filter(Gene %in% genes_highlight),
    size = 2,
    colour = "green",
  ) +
  
  # Resaltar genes
  geom_point(
    data = df_long %>% filter(type=="specific"),
    size = 2,
    colour = "blue", alpha = 0.4
  ) +
  
  # Etiquetas
  geom_text_repel(
    data = df_long %>% filter(Gene %in% genes_highlight),
    aes(label = Gene),
    size = 4,
    show.legend = FALSE,
    colour = "black"
  ) +
  
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  labs(
    title = "Relación entre Score y valor por métrica"
  )
g

#ggsave("~/tblab/yolanda/GLOWgenes/panelAPP/revision/degree_np_weight_plot.png", plot = g)

unique(df_long$metric)

orden_metric <- c(
  "functionGOBP",
  "coexpressionSTRING",
  "coessentAVANAwEXT",
  "fitnessRNAIEXT"
)

df_long_filtered <- df_long %>% 
  dplyr::filter(metric %in% orden_metric) %>%
  dplyr::filter(!(metric == "fitnessRNAIEXT" & value == 0)) %>%
  dplyr::filter(!(metric == "fitnessAVANAAEXT" & value == 0)) %>%
  dplyr::filter(value != 0) %>%
  dplyr::mutate(
    metric = factor(metric, levels = orden_metric)
  )

ggplot(df_long_filtered, aes(x = Score , y = value)) +
  geom_point(alpha = 0.3, colour = "red") +
  
  # Resaltar genes
  geom_point(
    data = df_long_filtered %>% filter(Gene %in% genes_highlight),
    size = 2,
    colour = "green", alpha = 0.8
  ) +
  
  # Resaltar genes específicos
  geom_point(
    data = df_long_filtered %>% filter(type == "specific"),
    size = 2,
    colour = "blue", alpha = 0.8
  ) +
  
  # Etiquetas
  geom_text_repel(
    data = df_long_filtered %>% filter(Gene %in% genes_highlight),
    aes(label = Gene),
    size = 4,
    show.legend = FALSE,
    colour = "black"
  ) +
  
  facet_wrap(~ metric, scales = "free") +
  theme_minimal(base_size = 14) +
  labs(
    title = NULL,
    x = "SGDS", 
    y = "Degree"
  )















