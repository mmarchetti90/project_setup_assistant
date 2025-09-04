#!/usr/bin/env python3

"""
This script checks samples counts distribution, and runs PCA and correlation analysis
"""

### ---------------------------------------- ###
### ANALYSIS TOOLS                           ###
### ---------------------------------------- ###

def cluster_data(data, n_neighbors=5):
    
    # Setting random seed (helps with clustering consistency)
    
    random.seed(42)
    
    # Computing kneighbors sparse matrix
    
    kneighbors_matrix = kneighbors_graph(X=data, n_neighbors=n_neighbors)
    
    # Creating edges list
    
    sources, targets = kneighbors_matrix.nonzero()
    edges = list(zip(sources, targets))
    
    # Building igraph object
    
    graph = igraph.Graph(directed=True)
    graph.add_vertices(kneighbors_matrix.shape[0])
    graph.add_edges(edges)
    
    # Converting graph to undirected
    
    graph.to_undirected(mode='collapse', combine_edges=None)
    
    # Clustering using Leiden algorithm
    
    clusters = graph.community_leiden(objective_function='modularity', weights=None, resolution_parameter=1.0, beta=0.01, initial_membership=None, n_iterations=2, node_weights=None).membership
    
    clusters = pd.DataFrame({'sample' : pca_data.index.values,
                             'condition' : clusters})
    
    return clusters

### ---------------------------------------- ###

def reduce_dimensions(cnts):
    
    # Fit PCA
    
    pca_model = PCA()
    pca_data = pca_model.fit_transform(cnts.T)
    
    pca_data = pd.DataFrame(pca_data, index=cnts.columns, columns=[f'PC{i+1}' for i in range(pca_data.shape[1])])
    
    # Extract explained variance
    
    explained_variance = pca_model.explained_variance_
    explained_variance = explained_variance * (100 / explained_variance.sum())
    explained_variance = explained_variance[:20]
    
    # Extract optimal components
    
    optimal_components, _ = kneedle(explained_variance)
    optimal_components += 1
    
    # Plot explained variance
    
    plt.figure(figsize=(7, 4))
    plt.plot([str(i+1) for i in range(explained_variance.shape[0])], explained_variance / 100, 'b', marker='o', markersize=5, linewidth=1)
    plt.vlines(optimal_components + 0.5, 0, max(explained_variance) / 100, linestyle='dashed', color='red', linewidth=1)
    plt.text(optimal_components + 1, max(explained_variance) / 110, 'Optimal PCA components', color='red')
    plt.xlabel('PC')
    plt.ylabel('Explained Variance (%)')
    plt.tight_layout()
    plt.savefig('pca_explained_variance.png', bbox_inches='tight', dpi=600)
    plt.close()
    
    # Fit UMAP
    
    umap_model = umap.UMAP(n_components=2, n_neighbors=5, random_state=42)

    umap_model.fit(pca_data)

    umap_data = umap_model.transform(pca_data)
    
    umap_data = pd.DataFrame(umap_data, index=pca_data.index, columns=['UMAP1', 'UMAP2'])
    
    # Export data
    
    pca_data.to_csv('pca_transform.tsv', sep='\t', index=True, header=True)
    
    umap_data.to_csv('umap_transform.tsv', sep='\t', index=True, header=True)
    
    pk.dump(pca_model, open('pca.pkl', 'wb'))
    
    pk.dump(umap_model, open('umap.pkl', 'wb'))
    
    return pca_model, pca_data, explained_variance, optimal_components, umap_model, umap_data

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
    
    plot_correlation_matrix(plot_data)

### ---------------------------------------- ###
### UTILS                                    ###
### ---------------------------------------- ###

def kneedle(vector, sort_vector=True):
    
    """
    Kneedle to find threshold cutoff.
    """
    
    if sort_vector:
        
        vector = np.sort(vector)[::-1]
    
    # Find gradient and intercept
    
    x0, x1 = 0, len(vector)
    y0, y1 = max(vector), min(vector)
    gradient = (y1 - y0) / (x1 - x0)
    intercept = y0
    
    # Compute difference vector
    
    difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(vector)]
    
    # Find max of difference_vector and define cutoff
    
    cutoff_index = difference_vector.index(max(difference_vector))
    cutoff_value = vector[cutoff_index]
    
    return cutoff_index, cutoff_value

### ---------------------------------------- ###

def normalize_counts(cnts):
    
    cnts = np.log1p(cnts * (1e6 / cnts.sum(axis=0)))
    
    return cnts

### ---------------------------------------- ###
### PLOTTING UTILS                           ###
### ---------------------------------------- ###

def plot_correlation_matrix(plot_data):
    
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
                            figsize=(len(unique_categories) * 1.5, 10))
        plt.xticks(rotation=90, fontweight='bold')
        plt.xlabel('Group', fontweight='bold')
        plt.ylabel('Frequency', fontweight='bold')
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        plt.savefig(out_name, bbox_inches='tight', dpi=600)
        plt.close()

### ---------------------------------------- ###

def plot_pca(pca_data, explained_variance, info=[], info_name='', out_plot='pca.png'):
    
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
        plot_data.loc[:, info_name] = [info.loc[info['sample'] == sample, 'condition'].values[0] if sample in info['sample'].values else 'NA' for sample in pca_data.index]
        
        # Plot data
        
        plt.figure(figsize=(5, 5))
        ax = sns.scatterplot(data=plot_data,
                             x='PC1',
                             y='PC2',
                             hue=info_name,
                             palette='tab10',
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

def plot_umap(umap_data, info=[], info_name='', out_plot='umap.png'):
    
    if len(info) == 0:
        
        # Prepare data for plotting
        
        plot_data = umap_data.loc[:, ['UMAP1', 'UMAP2']].copy()
        
        # Plot data
        
        plt.figure(figsize=(5, 5))
        ax = sns.scatterplot(data=plot_data,
                             x='UMAP1',
                             y='UMAP2',
                             s=50,
                             marker='o',
                             edgecolors='black',
                             linewidths=1)
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        plt.savefig(out_plot, bbox_inches='tight', dpi=600)
        plt.close()
    
    else:
        
        # Prepare data for plotting
        
        plot_data = umap_data.loc[:, ['UMAP1', 'UMAP2']].copy()
        plot_data.loc[:, info_name] = [info.loc[info['sample'] == sample, 'condition'].values[0] if sample in info['sample'].values else 'NA' for sample in pca_data.index]
        
        # Plot data
        
        plt.figure(figsize=(5, 5))
        ax = sns.scatterplot(data=plot_data,
                             x='UMAP1',
                             y='UMAP2',
                             hue=info_name,
                             palette='tab10',
                             s=50,
                             marker='o',
                             edgecolors='black',
                             linewidths=1)
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        plt.savefig(out_plot, bbox_inches='tight', dpi=600)
        plt.close()

### ------------------MAIN------------------ ###

import igraph
import numpy as np
import pandas as pd
import pickle as pk
import random
import seaborn as sns
import umap

from matplotlib import pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
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

    sample_info_files = []

### Number of detected genes

if '--min_counts' in argv:

    min_counts = int(argv[argv.index('--min_counts') + 1])

else:

    min_counts = 10

detected_genes = (counts > min_counts).sum(axis=0).reset_index(drop=False)

detected_genes.columns = ['sample', 'detected_genes']

detected_genes.loc[:, 'tot_sample_reads'] = counts.sum(axis=0).values

detected_genes.to_csv('detected_genes.tsv', sep='\t', index=False, header=True)

# Plot relation between tot_sample_reads, detected_genes

plot_data = np.log10(detected_genes[['tot_sample_reads', 'detected_genes']])

slope, intercept = np.polyfit(plot_data['tot_sample_reads'], plot_data['detected_genes'], 1)
trendline_x = plot_data['tot_sample_reads'].values
trendline_y_fun = np.poly1d((slope, intercept))
r, p = spearmanr(plot_data['tot_sample_reads'], plot_data['detected_genes'])

plot_data.plot(kind='scatter', x='tot_sample_reads', y='detected_genes', loglog=False, figsize=(5.1, 5))
plt.plot(trendline_x, trendline_y_fun(trendline_x), color='red', linestyle='--', label='Trendline')
plt.title(f'r = {r:.3f}\np = {p:.3e}', loc='left')
plt.xlabel('log$_{10}$(Sample reads)')
plt.ylabel('log$_{10}$(Detected genes)')
plt.savefig('detected_genes.png', bbox_inches='tight', dpi=600)
plt.close()

### Normalizing counts

print('Normalizing counts')

norm_counts = counts.loc[(counts > min_counts).sum(axis=1) >= 3,]

norm_counts = normalize_counts(norm_counts)

### Correlation

print('Running correlation analysis')

### PCA + UMAP

print('Reducing dimensions')

run_correlation(norm_counts)

pca_model, pca_data, explained_variance, optimal_components, umap_model, umap_data = reduce_dimensions(norm_counts)

print(f'Optimal PCA components: {optimal_components}')

plot_pca(pca_data, explained_variance, [], out_plot='pca.png')

plot_umap(umap_data, [], out_plot='umap.png')

### Cluster data based on PCA

print('Clustering data')

clusters = cluster_data(pca_data, n_neighbors=5)

clusters.to_csv('clusters.tsv', sep='\t', index=False)

sample_info_files.append('clusters.tsv')

### PCA/UMAP plots with sample info

for file in sample_info_files:

    category = file.split('/')[-1].replace('.tsv', '').replace('.txt', '.png')

    sample_info = pd.read_csv(file, sep='\t', index_col=None, header=0)

    out_name = 'pca_' + category + '.png'

    plot_pca(pca_data, explained_variance, sample_info, info_name=category, out_plot=out_name)
    
    out_name = 'umap_' + category + '.png'

    plot_umap(umap_data, sample_info, info_name=category, out_plot=out_name)

### Counts distribution

print('Checking counts distribution')

plot_counts_distribution(counts, sample_info_files)
