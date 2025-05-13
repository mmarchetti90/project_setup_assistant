#!/usr/bin/env python3

"""
This script checks samples counts distribution, and runs PCA and correlation analysis
"""

### ---------------------------------------- ###

def normalize_counts(cnts):
    
    cnts = np.log2(cnts * (1e6 / cnts.sum(axis=0)) + 1)
    
    return cnts

### ---------------------------------------- ###

def plot_counts_distribution(cnts, categories):
    
    # Plot counts distribution
    plot_data = pd.DataFrame({'Sample' : [col for col in cnts.columns for _ in range(cnts.shape[0])],
                              'Counts' : np.log10(np.concatenate([cnts.loc[:, col].values for col in cnts.columns]) + 1)})
    
    plt.figure(figsize=(cnts.shape[1], 10))
    ax = sns.boxplot(plot_data,
                     x='Sample',
                     y='Counts',
                     hue='Sample',
                     orient='v',
                     width=0.5,
                     whis=(0, 100),
                     dodge=False)
    plt.xticks(rotation=90, fontweight='bold')
    plt.xlabel(None)
    plt.ylabel('log10(Counts + 1)', fontweight='bold')
    #sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig('counts_distribution.png', bbox_inches='tight', dpi=600)
    plt.close()
    
    # Plot frequency of 0 counts in each sample
    plot_data = pd.DataFrame({"0" : [(cnts.loc[:, col].values == 0).sum() / cnts.shape[0] for col in cnts.columns],
                              "1+" : [(cnts.loc[:, col].values != 0).sum() / cnts.shape[0] for col in cnts.columns]},
                             index=cnts.columns.values)
    
    ax = plot_data.plot(kind='bar',
                        stacked=True,
                        color=['tomato', 'limegreen'],
                        edgecolor='black',
                        figsize=(cnts.shape[1], 10))
    plt.xticks(rotation=90, fontweight='bold')
    plt.xlabel(None)
    plt.ylabel('Frequency', fontweight='bold')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig('zero_frequency.png', bbox_inches='tight', dpi=600)
    plt.close()
    
    # Plot frequency of 0 counts in sample groups
    for file in categories:

        if file == []:
            
            continue

        out_name = 'zero_frequency_in_' + file.split('/')[-1].replace('.tsv', '.png').replace('.txt', '.png')

        cat = pd.read_csv(file, sep='\t', index_col=None, header=0)
        
        unique_categories = list(set(cat.condition.values))
        
        max_n = max([(cat.condition == c).sum() for c in unique_categories])
        
        plot_data = {i : [] for i in range(max_n + 1)}
        for c in unique_categories:
            
            c_samples = cat.loc[cat.condition == c, 'sample'].values
            
            zeros = (cnts.loc[:, cnts.columns.isin(c_samples)] == 0).sum(axis=1)
    
            for i in range(max_n + 1):
                
                plot_data[i].append((zeros == i).sum() / len(zeros))
        
        plot_data = pd.DataFrame(plot_data, index=unique_categories)
    
        ax = plot_data.plot(kind='bar',
                            stacked=True,
                            edgecolor='black',
                            figsize=(cnts.shape[1], 10))
        plt.xticks(rotation=90, fontweight='bold')
        plt.xlabel('Group', fontweight='bold')
        plt.ylabel('Frequency', fontweight='bold')
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        plt.savefig(out_name, bbox_inches='tight', dpi=600)
        plt.close()

### ---------------------------------------- ###

def run_pca(cnts):
    
    # Fit PCA
    
    pca = PCA()
    pca_data = pca.fit_transform(cnts.T)
    
    pca_data = pd.DataFrame(pca_data, index=cnts.columns, columns=[f'PC{i+1}' for i in range(pca_data.shape[1])])
    
    # Extract explained variance
    
    explained_variance = pca.explained_variance_
    explained_variance = explained_variance * (100 / explained_variance.sum())
    explained_variance = explained_variance[:20]
    
    # Plot explained variance
    
    plt.figure(figsize=(10, 4))
    plt.bar([str(i+1) for i in range(explained_variance.shape[0])],
            explained_variance,
            width=0.8,
            color='deepskyblue',
            edgecolor='black',
            linewidth=1)
    plt.xlabel('PC')
    plt.savefig('pca_explained_variance.png', bbox_inches='tight', dpi=600)
    plt.close()
    
    # Export data
    
    pca_data.to_csv('pca_transform.tsv', sep='\t', index=True, header=True)
    
    return pca_data, explained_variance
    
### ---------------------------------------- ###
    
def plot_pca(pca_data, explained_variance, info=[], out_plot='pca.png'):
    
    if len(info) == 0:
        
        # Prepare data for plotting
        plot_data = pca_data.loc[:, ['PC1', 'PC2']].copy()
        
        # Plot data
        plt.figure(figsize=(5, 5))
        ax = sns.scatterplot(data=plot_data,
                             x='PC1',
                             y='PC2',
                             s=50,
                             marker='o',
                             edgecolors='black',
                             linewidths=1)
        plt.xlabel(f'PC1 ({round(explained_variance[0], 2)}%)')
        plt.ylabel(f'PC2 ({round(explained_variance[1], 2)}%)')
        plt.savefig(out_plot, bbox_inches='tight', dpi=600)
        plt.close()
    
    else:
        
        # Prepare data for plotting
        plot_data = pca_data.loc[:, ['PC1', 'PC2']].copy()
        plot_data = plot_data.assign(Condition = [info.loc[info['sample'] == sample, 'condition'].values[0] if sample in info['sample'].values else 'NA' for sample in pca_data.index])
        
        # Plot data
        plt.figure(figsize=(5, 5))
        ax = sns.scatterplot(data=plot_data,
                             x='PC1',
                             y='PC2',
                             hue='Condition',
                             s=50,
                             marker='o',
                             edgecolors='black',
                             linewidths=1)
        plt.xlabel(f'PC1 ({round(explained_variance[0], 2)}%)')
        plt.ylabel(f'PC2 ({round(explained_variance[1], 2)}%)')
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        plt.savefig(out_plot, bbox_inches='tight', dpi=600)
        plt.close()

### ---------------------------------------- ###

def run_correlation(cnts):
    
    # Init results matrix
    plot_data = pd.DataFrame([[1. for _ in range(cnts.shape[1])] for _ in range(cnts.shape[1])],
                             index=cnts.columns,
                             columns=cnts.columns)
    
    # Compute correlation scores
    for a in range(cnts.shape[1]):
        
        for b in range(a, cnts.shape[1], 1):
            
            if a == b:
                
                continue
            
            #correlation, pval = pearsonr(cnts.iloc[:, a], cnts.iloc[:, b])
            correlation, pval = spearmanr(cnts.iloc[:, a], cnts.iloc[:, b])
            
            plot_data.iloc[a, b] = correlation
            plot_data.iloc[b, a] = correlation
            
            # Plot
            #plt.figure(figsize=(5, 5))
            #plt.plot(cnts.iloc[:, a],
            #         cnts.iloc[:, b],
            #         marker='.',
            #         markersize=5,
            #         lw=0)
            #plt.title(f'r = {round(correlation, 4)}\np-value = {round(pval, 4)}')
            #plt.xlabel(plot_data.columns[a])
            #plt.ylabel(plot_data.columns[b])
            #plt.savefig(f'{plot_data.columns[a]}_vs_{plot_data.columns[b]}.png', dpi=600)
            #plt.close()
    
    # Save data
    plot_data.to_csv('correlation_scores.tsv', sep='\t', index=True, header=True)

    # Plot
    plt.figure(figsize=(5, 5))
    sns.heatmap(plot_data,
                vmax=1.,
                cmap='crest',
                square=True,
                linecolor='black',
                linewidths=1)
    plt.tight_layout()
    plt.savefig('correlation.png', dpi=600)
    plt.close()
    
### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sys import argv

### Load data

counts_file = argv[argv.index('--counts_file') + 1]

counts = pd.read_csv(counts_file, sep='\t', index_col=None, header=0).iloc[:, 1:]

counts = counts.loc[:, np.sort(counts.columns.values)]

counts = counts.loc[(counts == 0).sum(axis=1) != counts.shape[1],]

### Load sample info

if '--sample_info_files' in argv:

    sample_info_files = argv[argv.index('--sample_info_files') + 1].split(',')

else:

    sample_info_files = [[]]

### Counts distribution

print('Checking counts distribution')

plot_counts_distribution(counts, sample_info_files)

### PCA

print('Running PCA')

norm_counts = counts.loc[(counts > 10).sum(axis=1) >= 3,]

norm_counts = normalize_counts(counts)

pca_data, explained_variance = run_pca(norm_counts)

plot_pca(pca_data, explained_variance, [], out_plot='pca.png')

### PCA plots

for file in sample_info_files:

    if file != []:

        out_name = 'pca_' + file.split('/')[-1].replace('.tsv', '.png').replace('.txt', '.png')

        sample_info = pd.read_csv(file, sep='\t', index_col=None, header=0)

        plot_pca(pca_data, explained_variance, sample_info, out_plot=out_name)

### Correlation

print('Running correlation analysis')

run_correlation(norm_counts)
