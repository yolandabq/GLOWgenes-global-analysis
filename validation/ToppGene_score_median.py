#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 16:31:41 2024

@author: yolanda
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



# Data.tsv is stored locally in the 
# same directory as of this python file
df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGenes_clustered_no_NA.tsv',sep = '\t') 
df.head()

df_np = df.iloc[:,:-1].set_index('SYMBOL').to_numpy()
df_np
panel_class = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/final_panels_list_2.txt',sep = '\t', header=None) 



def normalize_specificity_scores(specificity_scores):
    """Normalize the specificity scores to be between 0 and 1."""
    min_score = np.min(specificity_scores)
    max_score = np.max(specificity_scores)
    normalized_scores = (specificity_scores - min_score) / (max_score - min_score)
    return normalized_scores


masked_rank_matrix = np.ma.masked_equal(df_np, 0)
min_rankings = np.min(masked_rank_matrix, axis=1)
max_rankings = np.max(df_np, axis=0)

# Step 2: Calculate the median ranking for each element
#median_rankings = np.median(masked_rank_matrix, axis=1)
median_rankings = np.array([np.median(row.compressed()) for row in masked_rank_matrix])


# Step 3: Compute the specificity score (difference between median and minimum)
specificity_score = median_rankings - min_rankings


specificity_score_normalized = normalize_specificity_scores(specificity_score)


d_class={"Symbols": df["SYMBOL"], "Score": list(specificity_score_normalized)}

df_score_2=pd.DataFrame(data=d_class)
df_score_2["min_ranking"] = min_rankings
df_score_2["median_ranking"] = median_rankings

df_score_sorted_2=df_score_2.sort_values(by='Score', ascending=False)


plt.hist(list(df_score_sorted_2["Score"]), color='lightgreen', ec='black')

plt.scatter(median_rankings, min_rankings, s=5)
plt.show()


df_score_sorted_2.to_csv("/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/revision/ToppGene_score_median.tsv", index = False)






## nada, con la corrección esta queda mal. 
from collections import Counter

group_labels = np.array(panel_class[0])

# Step 4: Group-based adjustment
group_specific_scores = np.zeros(df_np.shape[0])

for i in range(df_np.shape[0]):
    # Calculate the number of top rankings in each group
    
    # Calculate the number of high rankings within each group
    group_counts = Counter()
    
    for j in range(df_np.shape[1]):
        if df_np[i, j] <= min_rankings[i]:  # Adjust threshold as needed
            group_counts[group_labels[j]] += 1
            print(group_counts)
    # Score adjustment: Favor concentration in fewer groups
    non_zero_groups = len(group_counts)
    group_specific_scores[i] = np.sum(np.array(list(group_counts.values())) ** 2) / (non_zero_groups if non_zero_groups > 0 else 1)


# Step 5: Adjust the specificity score with the group-specific score
adjusted_specificity_score = specificity_score * group_specific_scores

adjusted_specificity_score_normalized = normalize_specificity_scores(adjusted_specificity_score)


d_class={"Symbols": df["SYMBOL"], "Score": list(adjusted_specificity_score_normalized)}

df_adjusted_score_2=pd.DataFrame(data=d_class)

df_adjusted_score_sorted_2=df_adjusted_score_2.sort_values(by='Score', ascending=False)

df_adjusted_score_sorted_2[df_adjusted_score_sorted_2["Symbols"] == "MAJIN"]

plt.hist(list(df_adjusted_score_sorted_2["Score"]), color='lightgreen', ec='black')


import matplotlib.pyplot as plt
plt.scatter(list(median_rankings),list(df_score_2["Score"]) )



