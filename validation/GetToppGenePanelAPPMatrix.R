#!/usr/bin/env Rscript

library("dplyr")

input_toppgene_dir <- c("/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/graci_output_renamed/", 
               "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/pablo_output/", 
               "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/cris_output/", 
               "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/yoli_output/")

#glowed_panels <- list.files(paste(input_dir,"/results_GA/", sep = ""))
#"/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list.txt"
glowed_panels_list <-  read.delim("/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list.txt", header = F, stringsAsFactors = F, quote = "", check.names=F)
glowed_panels <- sub("$","_GA",glowed_panels_list$V2)

## obtenemos los rankings de ToppGene

toppgene_list <- c()

for (f in input_toppgene_dir) {
  toppgene_results <- list.files(
    path = f,
    pattern = "*_GA.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  toppgene_list <- c(toppgene_list, toppgene_results)
}

### vamos a obtener también los paneles de genes originales (input de ToppGene)

PanelApp_panels_list <- c()

input_panels_dir <- c("/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/graci/", 
                        "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/yoli/", 
                        "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/cris/", 
                        "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/panels_for_toppgene/pablo/")

for (f in input_panels_dir) {
  panels <- list.files(
    path = f,
    pattern = "*_GA.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  PanelApp_panels_list <- c(PanelApp_panels_list, panels)
}

df_panels <- data.frame(
  panel = tools::file_path_sans_ext(basename(PanelApp_panels_list)),
  path  = PanelApp_panels_list,
  stringsAsFactors = FALSE
)

#############

df_glowmatrix=data.frame(SYMBOL = character()) 
n_panel <- 1

for (g_panel in toppgene_list){
  panel_name <- tools::file_path_sans_ext(basename(g_panel))
  print(panel_name)
  
  n_panel <- n_panel +1
  #print(paste(input_dir,"/results/", g_panel, sep = ""))
  glowgenes = read.delim(g_panel, header = T, stringsAsFactors = F, check.names=F, sep = ",")
  
  ranking = glowgenes[c("GeneSymbol", "Rank")]
  colnames(ranking) <- c("SYMBOL", "Rank")
  
  panel=read.delim(df_panels[df_panels$panel == panel_name, "path"], header = T, stringsAsFactors = F, quote = "", check.names=F)
  colnames(panel) <- "SYMBOL"
  panel$Rank=0

  ranking <- rbind(ranking, panel)
  
  colnames(ranking)[2] <- panel_name
  
  rownames(ranking) <- ranking$SYMBOL 
  
  print(c("Nº genes antes del merge:", nrow(df_glowmatrix)))
  df_glowmatrix <- dplyr::full_join(
    df_glowmatrix,ranking[c( "SYMBOL",panel_name)],by="SYMBOL"
  )
  print(c("Nº genes después del merge:", nrow(df_glowmatrix)))
  print("---------------------------") 
}


write.table(df_glowmatrix, "/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGenes.tsv",sep="\t", row.names=FALSE, quote=FALSE)


#panel=read.delim("/mnt/tblab/yolanda/GLOWgenes/Panel_DHR.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) # cuidado! comprobar si el panel tiene cabecera o no
#glowgenes = read.delim("/mnt/tblab/yolanda/GLOWgenes/Retinal_disorders_GA/GLOWgenes_prioritization_Random.txt", header = F, stringsAsFactors = F, quote = "", check.names=F)
#colnames(glowgenes) = c("SYMBOL", "score_GA", "GLOWgenes_GA")
#colnames(panel) = "genes"
