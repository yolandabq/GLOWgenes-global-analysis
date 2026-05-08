#!/usr/bin/env Rscript

library("dplyr")

input_dir <- "/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis"
#glowed_panels <- list.files(paste(input_dir,"/results_GA/", sep = ""))
#"/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list.txt"
glowed_panels_list <-  read.delim("/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list.txt", header = F, stringsAsFactors = F, quote = "", check.names=F)
glowed_panels <- sub("$","_GA",glowed_panels_list$V2)

df_glowmatrix=data.frame(SYMBOL = character()) 
n_panel <- 1
for (g_panel in glowed_panels){
  print(c(g_panel, " ,panel ", n_panel))
  n_panel <- n_panel +1
  #print(paste(input_dir,"/results/", g_panel, sep = ""))
  glowgenes = read.delim(paste(input_dir,"/results_GA/", g_panel, "/GLOWgenes_prioritization_Random.txt", sep = ""), header = F, stringsAsFactors = F, quote = "", check.names=F)
  colnames(glowgenes) = c("SYMBOL", "score", "GLOWgenes")
  
  panel=read.delim(paste(input_dir,"/panels/", g_panel, ".csv", sep=""), header = T, stringsAsFactors = F, quote = "", check.names=F)
  colnames(panel) <- "SYMBOL"
  panel$score=0
  panel$GLOWgenes =0
  
  glowgenes <- rbind(glowgenes, panel)
  
  colnames(glowgenes)[3] <- g_panel
  
  rownames(glowgenes) <- glowgenes$SYMBOL 

  print(c("Nº genes antes del merge:", nrow(df_glowmatrix)))
  df_glowmatrix <- dplyr::full_join(
    df_glowmatrix,glowgenes[c( "SYMBOL",g_panel)],by="SYMBOL"
  )
  print(c("Nº genes después del merge:", nrow(df_glowmatrix)))
  print("---------------------------") 
}


write.table(df_glowmatrix, "/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/glowmatrix.tsv",sep="\t", row.names=FALSE, quote=FALSE)


#panel=read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/Panel_DHR.txt", header = F, stringsAsFactors = F, quote = "", check.names=F) # cuidado! comprobar si el panel tiene cabecera o no
#glowgenes = read.delim("/home/yolanda/tblab/yolanda/GLOWgenes/Retinal_disorders_GA/GLOWgenes_prioritization_Random.txt", header = F, stringsAsFactors = F, quote = "", check.names=F)
#colnames(glowgenes) = c("SYMBOL", "score_GA", "GLOWgenes_GA")
#colnames(panel) = "genes"
