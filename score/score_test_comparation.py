#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 17:54:42 2024

@author: yolanda
"""

import pandas as pd
import scipy.stats as stats

# Step 1: Load your data
# Assuming the main table with symbols and scores is already in a CSV or similar format.
# You should replace 'main_table.csv', 'class1_symbols.csv', 'class2_symbols.csv' with actual file paths
main_df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/score/score_median.tsv')  # Your main table with symbols and scores
panel_specific_df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_panel_specific.txt')  # Your table with Class 1 symbols
#class_specific_df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_class_specific_05.txt')  # Your table with Class 2 symbols

#class_specific_df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_class_specific.txt')  # Your table with Class 2 symbols
#class_specific_df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_class_specific_top1_bottomno0.txt')  # Your table with Class 2 symbols
class_specific_df = pd.read_csv('/home/yolanda/tblab/yolanda/GLOWgenes/panelAPP/analysis/genes_specific/genes_class_specific_laxo.txt')  # Your table with Class 2 symbols


# Step 2: Separate the elements into groups (Class1, Class2, Others)
panel_specific_symbols = set(panel_specific_df['SYMBOL'])  # Replace 'Symbols' with actual column name for class1
class_specific_symbols = set(class_specific_df['SYMBOL'])  # Replace 'Symbols' with actual column name for class2

# Create new columns to classify the symbols in the main table
def classify_symbol(symbol):
    if symbol in panel_specific_symbols:
        return 'Panel_specific'
    elif symbol in class_specific_symbols:
        return 'Class_specific'
    else:
        return 'Others'

main_df['Class'] = main_df['Symbols'].apply(classify_symbol)

# Step 3: Split the data based on class
panel_specific_scores = main_df[main_df['Class'] == 'Panel_specific']['Score']
class_specific_scores = main_df[main_df['Class'] == 'Class_specific']['Score']
others_scores = main_df[main_df['Class'] == 'Others']['Score']

#panel_specific_scores=panel_specific_scores[panel_specific_scores!="--"]
#class_specific_scores=class_specific_scores[class_specific_scores!="--"]
#others_scores=others_scores[others_scores!="--"]

# Step 4: Perform statistical analysis
# Here we perform ANOVA to see if there are differences in means between groups
anova_result = stats.f_oneway(panel_specific_scores, class_specific_scores, others_scores)
print(f"ANOVA Result: F-statistic = {anova_result.statistic}, p-value = {anova_result.pvalue}")

# If you want to use a non-parametric test (e.g., Kruskal-Wallis), you can use this instead
kruskal_result = stats.kruskal(panel_specific_scores, class_specific_scores, others_scores)
print(f"Kruskal-Wallis Result: H-statistic = {kruskal_result.statistic}, p-value = {kruskal_result.pvalue}")



import seaborn as sns
import matplotlib.pyplot as plt

# Create the boxplot to visualize the distributions of scores across the classes
plt.figure(figsize=(8,6))
sns.boxplot(x='Class', y='Score', data=main_df)
plt.title('Score Distributions by Class')
plt.xlabel('Class')
plt.ylabel('Score')
plt.show()



plt.hist(list(panel_specific_scores), color='lightgreen', ec='black', bins=20)
plt.show()

plt.hist(list(class_specific_scores), color='lightgreen', ec='black', bins=20)
plt.show()


plt.hist(list(others_scores), color='lightgreen', ec='black', bins=50)
plt.show()


#---------


from itertools import combinations
from scipy.stats import mannwhitneyu

# Sample 400 random values from "Others" group
others_sample = main_df[main_df['Class'] == 'Others'].sample(n=400, random_state=42)

# Create the new groups for analysis

panel_specific_scores = main_df[main_df['Class'] == 'Panel_specific']['Score']
class_specific_scores = main_df[main_df['Class'] == 'Class_specific']['Score']
others_scores = others_sample['Score']  # Use the random sample for "Others"



# Update the 'main_df' with the 400 random samples for "Others"
updated_df = pd.concat([main_df[main_df['Class'] != 'Others'], others_sample])


#panel_specific_scores=panel_specific_scores[panel_specific_scores!="--"]
#class_specific_scores=class_specific_scores[class_specific_scores!="--"]
#others_scores=others_scores[others_scores!="--"]

# Perform pairwise Mann-Whitney U tests with the random "Others" subset
group_scores = {
    'Panel_specific': panel_specific_scores,
    'Class_specific': class_specific_scores,
    'Others': others_scores  # Use the random sample for "Others"
}



for (group1, group2) in combinations(group_scores.keys(), 2):
    u_stat, p_value = mannwhitneyu(group_scores[group1], group_scores[group2], alternative='two-sided')
    print(f"Comparison between {group1} and {group2}: U-statistic = {u_stat}, p-value = {p_value}")


# Combine the class data with the random sample for "Others"
updated_df = pd.concat([main_df[main_df['Class'] != 'Others'], others_sample])

# Create the boxplot
plt.figure(figsize=(8,6))
sns.boxplot(x='Class', y='Score', data=updated_df)
plt.title('Score Distributions by Class (with Random Sample of Others)')
plt.xlabel('Class')
plt.ylabel('Score')
plt.show()



### cuantos más bottons, más score. 

class_specific = class_specific_df.merge(main_df, how="left", left_on='SYMBOL', right_on='Symbols')


plt.figure(figsize=(8,6))
sns.boxplot(x='n_bottom', y='Score', data=class_specific)
plt.title('Score Distributions by Class number bottom')
plt.xlabel('N_bottom')
plt.ylabel('Score')
plt.show()


### iteramos 10 veces el test y el boxplot con una muestra random de others. 

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import numpy as np
from itertools import combinations
from scipy.stats import mannwhitneyu
# Assume 'main_df' has already been defined, and contains the 'Symbols', 'Score', and 'Class' columns

# Get the fixed class1 and class2 scores
# Create the new groups for analysis

panel_specific_scores = main_df[main_df['Class'] == 'Panel_specific']['Score']
class_specific_scores = main_df[main_df['Class'] == 'Class_specific']['Score']


# Define the groups for the Mann-Whitney U test
group_labels = ['Panel_specific', 'Class_specific', 'Others']
# Store the results of pairwise tests across 10 iterations
pairwise_results = []

# Perform Kruskal-Wallis test and boxplot for 10 iterations
kruskal_results = []

for i in range(10):
    # Randomly sample 400 elements from "Others"
    others_sample = main_df[main_df['Class'] == 'Others'].sample(n=400, random_state=i)
    others_scores = others_sample['Score']  # Scores of the random sample from "Others"
    
    # Perform Kruskal-Wallis test
    kruskal_result = kruskal(panel_specific_scores, class_specific_scores, others_scores)
    kruskal_results.append({
        'Iteration': i+1,
        'H-statistic': kruskal_result.statistic,
        'p-value': kruskal_result.pvalue
    })
    
    # Prepare data for the boxplot
    updated_df = pd.concat([
        main_df[main_df['Class'] != 'Others'],  # All data except "Others"
        others_sample  # Add the 400 random "Others"
    ])
    
    # Generate the boxplot for this iteration
    plt.figure(figsize=(8,6))
    sns.boxplot(x='Class', y='Score', data=updated_df)
    plt.title(f'Score Distributions by Class (Iteration {i+1})')
    plt.xlabel('Class')
    plt.ylabel('Score')
    plt.show()
    
    
    # Perform pairwise Mann-Whitney U tests for each pair of groups
    
    group_scores = {
        'Panel_specific': panel_specific_scores,
        'Class_specific': class_specific_scores,
        'Others': others_scores  # Use the random sample for "Others"
    }

    
    for (group1, group2) in combinations(group_labels, 2):
        u_stat, p_value = mannwhitneyu(group_scores[group1], group_scores[group2], alternative='two-sided')
        pairwise_results.append({
            'Iteration': i + 1,
            'Group1': group1,
            'Group2': group2,
            'U-statistic': u_stat,
            'p-value': p_value
        })

# Convert Kruskal-Wallis results to a DataFrame
kruskal_df = pd.DataFrame(kruskal_results)

# Display Kruskal-Wallis results
print(kruskal_df)

# Convert the results into a DataFrame for easy viewing
pairwise_df = pd.DataFrame(pairwise_results)

# Display the pairwise Mann-Whitney U test results
print(pairwise_df)


plt.figure(figsize=(8, 6))
sns.violinplot(x='Class', y='Score', data=updated_df)
plt.title(f'Score Distributions by Class (Violin Plot)')
plt.xlabel('Class')
plt.ylabel('Score')
plt.show()



plt.figure(figsize=(8, 6))
sns.stripplot(x='Class', y='Score', data=updated_df, jitter=True)
plt.title(f'Score Distributions by Class (Strip Plot)')
plt.xlabel('Class')
plt.ylabel('Score')
plt.show()





sns.pairplot(main_df, hue='Class')
plt.show()









