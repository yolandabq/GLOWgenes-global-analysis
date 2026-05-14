#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:05:21 2026

@author: yolanda
"""

import pandas as pd
import os
import networkx as nx
import numpy as np

networks_path = "/mnt/tblab/yolanda/GLOWgenes/berta/networks/pSNOW_databases_new_nets/"

## cargamos los genes específicos

disease_specific_genes = pd.read_csv("/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_panel_specific_panels.txt", sep= ",")



## cargamos la tabla de stats tops y cambiamos el nombre de los clusters
stat_tops = pd.read_csv("/mnt/tblab/yolanda/GLOWgenes/panelAPP/analysis/data_exploration/stat_tops_vfinal.tsv", sep = "\t" )

mapping = {
    4: 1,
    3: 2,
    0: 3,
    2: 4,
    1: 5
}

stat_tops["cluster"] = stat_tops["cluster"].map(mapping).fillna(stat_tops["cluster"])


pd_networks_names = pd.read_csv("/mnt/tblab/yolanda/GLOWgenes/berta/networks/networks_paths.txt", sep="\t", header = None)
network_names = [os.path.basename(x) for x in pd_networks_names[0]]


#### ahora ordenamos por porc_tops y nos quedamos con el mismo número de genes específicos y generales

stat_tops = stat_tops.sort_values(by="Porc_panels_top", ascending=False)
general_genes = stat_tops[0:disease_specific_genes.shape[0]]
general = general_genes.rename(columns={"SYMBOL": "Gene"})


specific = disease_specific_genes.copy()
general["type"] = "general"
specific["type"] = "specific"
general = general[["Gene", "type"]]
specific = specific[["Gene", "type"]]
merged = pd.concat([general, specific], ignore_index=True) # genes generales y específicos

####### calculamos el betweenes en esos genes

# copia para no machacar
final_df = merged.copy()

genes = merged["Gene"].str.upper().unique()

for net_name in network_names:
#for net_name in ["functionGOBP_HGNCnets.txt", "phenotypeHPO_HGNCnets.txt", "coocurrenceSTRING_HGNCnets.txt", "fitnessAVANAAEXT_HGNCnets.txt", "coexpressionCOXPRESdbEXT_HGNCnets.txt"]:
#for net_name in ["functionGOBP_HGNCnets.txt", "phenotypeHPO_HGNCnets.txt"]:

    print(f"Procesando: {net_name}")
    
    # --- 1. Cargar red ---
    network = pd.read_csv(
        os.path.join(networks_path, net_name),
        sep=" ",
        header=None
    )
    
    network.columns = ["gene1", "gene2", "weight_raw", "extra"]
    
    # limpiar
    network["weight"] = network["extra"].str.extract(r"([0-9.]+)").astype(float)
    network = network[["gene1", "gene2", "weight"]]
    
    # normalizar nombres
    network["gene1"] = network["gene1"].str.upper()
    network["gene2"] = network["gene2"].str.upper()
    
    ## transformar el peso a una distancia inversa antes de calcular la betweenness
    network["distance"] = 1 / network["weight"]
    
    # --- 2. Grafo ---
    G = nx.from_pandas_edgelist(
        network,
        source="gene1",
        target="gene2",
        edge_attr=["weight", "distance"]
    )
    
    # --- 3. Filtrar genes presentes ---
    genes_in_graph = [g for g in genes if g in G]
    
    print(f"  Genes en red: {len(genes_in_graph)} / {len(genes)}")
    
    if len(genes_in_graph) == 0:
        continue  # evitar errores
    
    # --- 4. Degree ---
    #deg = dict(G.degree(weight="weight"))
    deg = dict(G.degree())
    
    
    deg_df = pd.DataFrame({
        "Gene": stat_tops["SYMBOL"],
        "degree": [deg.get(g, np.nan) for g in stat_tops["SYMBOL"]]
    })
    

    # --- 5. Nombre de columna (limpio) ---
    col_name = net_name.replace("_HGNCnets.txt", "")
    
    deg_df = deg_df.rename(columns={"degree": col_name})
    
    # --- 6. Merge ---
    final_df = final_df.merge(
        deg_df,
        on="Gene",
        how="outer"
    )
    
###########################


#final_df.to_csv("/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/bewteenness_specific_general_2_k50.csv", index=False) 
final_df.to_csv("/mnt/tblab/yolanda/GLOWgenes/panelAPP/revision/degree_no_weight.csv", index=False) 



