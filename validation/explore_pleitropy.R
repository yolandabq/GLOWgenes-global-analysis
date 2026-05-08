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


plot(go_bp_counts$n_GO_BP,
     go_bp_counts$Score,
     xlab = "Nº GOs",
     ylab = "SGDS",
     pch = 16)

cor(go_bp_counts$n_GO_BP, go_bp_counts$Score)

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
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

######################

sscore_glow <- read.table("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv", header = TRUE, sep = ",")
betweenness <- read.table("~/tblab/yolanda/GLOWgenes/panelAPP/revision/bewteenness_specific_general.csv", header = TRUE, sep = ",")

df_long <- betweenness %>%
  pivot_longer(
    cols = -c(Gene, type),
    names_to = "metric",
    values_to = "value"
  )

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

ggplot(df_long, aes(x = type, y = value_log, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal()

df_long %>% dplyr::filter(metric %in% c("functionGOBP", "textminingSTRING", "coessentAVANAwEXT" ,"fitnessAVANAAEXT")) %>%
  ggplot( aes(x = type, y = value_log, fill = type)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free") +
  labs(
    title = "Comparación General vs Specific por métrica",
    x = "GLOWgenes network",
    y = "Betweenness (log)"
  )+
  theme_minimal(base_size = 15)

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




