#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 17:14:42 2024

@author: yolanda
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 


# sacamos el % de genes de un panel que están incluidos en otro. 


def panels_inclusion(panel1, panel2): 
    common_elements = set(panel1).intersection(set(panel2))
    num_common_elements = len(common_elements)
     
    num_panel1_elements = len(set(panel1))
     
    # Calculate the percentage similarity
    percentage_inclusion = (num_common_elements / num_panel1_elements) * 100
    
    return percentage_inclusion



gene_list = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list_2.txt',sep = '\t', header=None) 

gene_list.columns = ["Class","Panel"]

gene_list


### descripción de los genes de cada panel

panel_genes = dict()

for panel in gene_list.Panel:
    #print(panel)
    #print(str("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/panels/" + panel + "_GA.csv"))
    
    genes = pd.read_csv(str("~/tblab/yolanda/GLOWgenes/panelAPP/analysis/panels/" + panel + "_GA.csv"),sep = '\t') 
    panel_genes[panel] = genes["gene_data.hgnc_symbol"]


# descripción de los paneles

n_genes_panel=dict()
n_genes=[]
for p in panel_genes.keys():
    print(p + " -> " + str(len(panel_genes[p])))
    n_genes.append(len(panel_genes[p]))

n_genes_panel = {'panels' : panel_genes.keys(), 'n_genes' : n_genes}

df_n_panel_genes = pd.DataFrame.from_dict(n_genes_panel)
df_n_panel_genes["Class"] = gene_list["Class"]


plt.figure(figsize=(8, 6))
g = sns.boxplot(y='Class', x='n_genes', data=df_n_panel_genes)
#g.set(xticklabels=[])
plt.title(f'Number of genes by Class')
plt.xlabel('Class')
plt.ylabel('Number of genes')
plt.show()


all_genes_panel = [ list(panel_genes[x]) for x in panel_genes.keys()]
all_genes = pd.DataFrame(list(set(sum(all_genes_panel, []))))
all_genes = all_genes.rename(columns = {0:"genes"})

all_genes.to_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/gene_sets/all_genes_panelapp.tsv', sep="\t", index=False)

######

inclusion_matrix = np.empty([len(panel_genes.keys()), len(panel_genes.keys())])
inclusion_matrix[0:,]

r = 0
for panel1 in gene_list.Panel:
    c = 0
    for panel2 in gene_list.Panel:
        print(panel1 + " - " + panel2)
        inclusion_matrix[r,c] = panels_inclusion(list(panel_genes[panel1]),list(panel_genes[panel2]))
        c = c+1
    r = r+1
        
# plotting the heatmap 
hm = sns.heatmap(data=inclusion_matrix) 
plt.show()

df_inclusion = pd.DataFrame(inclusion_matrix, columns=gene_list.Panel, index=gene_list.Panel)


############


for category in set(gene_list.Class):
    print(gene_list[gene_list["Class"] == category]["Panel"])
    print("--")
    panels_c = gene_list[gene_list["Class"] == category]["Panel"]
    

hm = sns.heatmap(data=inclusion_matrix) 
plt.show()



set(gene_list.Class)
panels_c = gene_list[gene_list["Class"] == "Endocrine disorders"]["Panel"]
p=df_inclusion.loc[panels_c, panels_c][df_inclusion.loc[panels_c, panels_c] == 100]
p=df_inclusion.loc[panels_c, panels_c][df_inclusion.loc[panels_c, panels_c] > 90]
p=df_inclusion.loc[panels_c, panels_c][df_inclusion.loc[panels_c, panels_c] > 50]

hm = sns.heatmap(data=df_inclusion.loc[panels_c, panels_c]) 
plt.show()


import plotly.express as px
import os

# Ensure the output directory exists
output_dir = "heatmaps"
os.makedirs(output_dir, exist_ok=True)

# Loop through the DataFrames and create heatmaps
for panel_class  in list(set(gene_list.Class)):
    panels_c = gene_list[gene_list["Class"] == panel_class]["Panel"]
    
    print(panels_c)
    fig = px.imshow(df_inclusion.loc[panels_c, panels_c],
                    x=df_inclusion.loc[panels_c, panels_c].columns, y=df_inclusion.loc[panels_c, panels_c].index,
                    title=f"Heatmap {panel_class}"
                   )

    # Save each figure as an HTML file
    fig.write_html(os.path.join(output_dir, f"heatmap_{panel_class}.html"))

    # Alternatively, save each figure as a PNG file
    fig.write_image(os.path.join(output_dir, f"heatmap_{panel_class}.png"))

    # If you want to display them in the loop
    # fig.show(renderer='browser')


fig = px.imshow(df_inclusion.loc[panels_c, panels_c], x=df_inclusion.loc[panels_c, panels_c].columns, y=df_inclusion.loc[panels_c, panels_c].index)
fig.update_xaxes(side="top")
fig.show()


# vamos a filtrar, cuando uno está al 100% en otro, quitamos el grande. 

remove_superpanel=[]
remove_superpanel_dict = {}

for panel_class  in list(set(gene_list.Class)):
    panels_c = gene_list[gene_list["Class"] == panel_class]["Panel"]
    panels_to_remove_class = []
    for panel_c_i in panels_c:
        for panel_c_c in panels_c:
            if df_inclusion.loc[panel_c_i, panel_c_c] == 100 and panel_c_i != panel_c_c: 
                remove_superpanel.append(panel_c_c)
                panels_to_remove_class.append(panel_c_c)
    remove_superpanel_dict[panel_class] = list(set(panels_to_remove_class))
    

unique_panels_remove = pd.DataFrame(list(set(remove_superpanel)))

unique_panels_remove.to_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/filter_panels/panels_including_others.txt', sep="\t", index=False)



















