#!/usr/bin/env python3

"""
This script reads in Kraken2/Bracken reports and merges them at the desired level
"""

### ---------------------------------------- ###

def parse_args():

    # Kraken2 files
    k2_reports = argv[argv.index('--k2_reports') + 1].split(',')
    
    # Sample IDs
    reports_ids = argv[argv.index('--reports_ids') + 1].split(',')

    # Taxonomic level
    # (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies
    tax_lvl = argv[argv.index('--level') + 1]
    
    # Contaminant species (species name or taxid)
    if '--contaminants' in argv:

        contaminants = argv[argv.index('--contaminants') + 1].split(',')
        contaminants = [c.replace('_', ' ') for c in contaminants]

    else:

        contaminants = []

    # Sample info table
    if '--sample_info' in argv:

        sample_info_path = argv[argv.index('--sample_info') + 1]

    else:

        sample_info_path = ''
    
    # PCA model path
    if '--pca_model' in argv:

        pca_model_path = argv[argv.index('--pca_model') + 1]

    else:

        pca_model_path = ''
    
    # UMAP model path
    if '--umap_model' in argv:

        umap_model_path = argv[argv.index('--umap_model') + 1]

    else:

        umap_model_path = ''
    
    # Rescale data toggle
    if '--rescale' in argv:
        
        rescale_toggle = True
    
    else:
        
        rescale_toggle = False

    # Use raw counts or percent for downstream analyses?
    # Usually best to use percent
    if '--use_raw' in argv:

        use_raw = True

    else:

        use_raw = False

    return k2_reports, reports_ids, tax_lvl, contaminants, sample_info_path, pca_model_path, umap_model_path, rescale_toggle, use_raw

### ---------------------------------------- ###

def parse_k2_report(sample_name, path, desired_taxonomy, contaminants=[]):
    
    k2_data_raw, k2_data_percent = {}, {}
    
    for line in open(path, 'r').readlines():
        
        # Trim line
        line = line.strip()
        
        if not len(line):
            
            continue
        
        line = line.replace('\n', '')
        
        # Extract info
        if len(line.split('\t')) == 8:
            
            percent_fragments_clade, n_fragments_clade, n_fragments_direct, taxon_minimizers, taxon_unique_minimizers, rank, taxid, name = line.split('\t')
        
        else:
            
            percent_fragments_clade, n_fragments_clade, n_fragments_direct, rank, taxid, name = line.split('\t')

        # Trim percent_fragments_clade
        percent_fragments_clade = percent_fragments_clade.strip()
        percent_fragments_clade = float(percent_fragments_clade)
        
        # Trim n_fragments_clade
        n_fragments_clade = n_fragments_clade.strip()
        n_fragments_clade = int(n_fragments_clade)

        # Trim name
        while name.startswith(' '):
            
            name = name[1:]
        
        # Skip if undesired
        root_n, unclassified_n = 0, 0
        if rank == 'R':
            
            root_n = percent_fragments_clade
        
        elif rank == 'U':
            
            unclassified_n = percent_fragments_clade
        
        elif rank not in desired_taxonomy or n_fragments_clade == 0 or taxid in contaminants or name in contaminants:
            
            continue
        
        else:
            
            k2_data_raw[name] = n_fragments_clade
            k2_data_percent[name] = percent_fragments_clade
        
    k2_data_raw = pd.Series(k2_data_raw, name=sample_name)
    k2_data_percent = pd.Series(k2_data_percent, name=sample_name)
    
    return k2_data_raw, k2_data_percent, root_n, unclassified_n

### ---------------------------------------- ###

def scale_features(raw_dt):
    
    # Copying normalized_counts matrix
    scaled_features = raw_dt.copy()
    
    # Scaling features
    features_mean = raw_dt.iloc[:, 1:].mean(axis = 1).to_numpy()
    features_std = raw_dt.iloc[:, 1:].std(axis = 1).to_numpy()
    scaled_features.iloc[:, 1:] = scaled_features.iloc[:, 1:].subtract(features_mean, axis="index").div(features_std, axis="index")
    
    # Deal with NAs
    scaled_features = scaled_features.fillna(0)
    
    return scaled_features

### ---------------------------------------- ###

def reduce_dimensions(dt, pca_model_path='', umap_model_path='', pca_components=50, neighbors=30):
    
    if pca_components > dt.shape[1]:

        pca_components = dt.shape[1]
    
    # Scale data
    dt_scaled = scale_features(dt)
    
    # PCA transform
    pca_input = dt_scaled.loc[:, ~ dt_scaled.columns.isin(['clade'])].T
    
    try:
        
        # Transform data
        pca_model = pk.load(open(pca_model_path, "rb"))
        pca_dt = pca_model.transform(pca_input)
        
        # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
        x0, x1 = 0, min(len(pca_model.explained_variance_), pca_components - 1)
        y0, y1 = pca_model.explained_variance_[x0], pca_model.explained_variance_[x1]
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(pca_model.explained_variance_[:pca_components])]
        optimal_pca_components = difference_vector.index(max(difference_vector)) + 1
    
    except:
        
        # Transform data
        pca_model = PCA(pca_components)
        pca_dt = pca_model.fit_transform(pca_input)
    
        # Selecting optimal number of PCs using the elbow method (simplified Kneedle)
        x0, x1 = 0, min(len(pca_model.explained_variance_), pca_components - 1)
        y0, y1 = pca_model.explained_variance_[x0], pca_model.explained_variance_[x1]
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(pca_model.explained_variance_[:pca_components])]
        optimal_pca_components = difference_vector.index(max(difference_vector)) + 1
    
        # Elbow plot of explained variance
        plt.figure(figsize=(5, 3))
        plt.plot(range(x0 + 1, x1 + 2), pca_model.explained_variance_[:50] / 100, 'b', marker='o', markersize=5, linewidth=1)
        plt.plot([optimal_pca_components, optimal_pca_components], [y0 / 100, y1 / 100], linestyle='dashed', color='red', linewidth=1)
        plt.xlabel('PC')
        plt.ylabel('Explained Variance (%)')
        plt.tight_layout()
        plt.savefig("PCA_ExplainedVariance.png", dpi=300)
        plt.close()
        
        # Save model
        pk.dump(pca_model, open("pca.pkl", "wb"))
    
    # UMAP transfrom
    try:
        
        # Transform data
        umap_model = pk.load(open(umap_model_path, "rb"))
        umap_dt = umap_model.transform(pca_dt)
    
    except:
        
        # Transform data
        umap_model = umap.UMAP(n_components=2, n_neighbors=neighbors, random_state=42)
        umap_dt = umap_model.fit_transform(pca_dt)
    
        # Save model
        pk.dump(umap_model, open("umap.pkl", "wb"))
    
    pca_dt = pd.DataFrame(data = pca_dt[:, :optimal_pca_components], index = dt.columns[1:], columns = [f'PC_{i+1}' for i in range(optimal_pca_components)])
    umap_dt = pd.DataFrame(data = umap_dt, index = dt.columns[1:], columns = ['UMAP_1', 'UMAP_2'])
    
    return pca_dt, umap_dt

### ---------------------------------------- ###

def cluster_samples(dt, n_neighbors=10):
    
    # Setting random seed (helps with clustering consistency)
    random.seed(42)
    
    # Computing kneighbors sparse matrix
    kneighbors_matrix = kneighbors_graph(X=dt.loc[:, dt.columns != 'sample'], n_neighbors=n_neighbors)
    
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
    
    # Convert to pandas data frame
    clusters = pd.DataFrame({'sample' : dt['sample'],
                             'cluster' : clusters})
    
    return clusters

### ---------------------------------------- ###

def plot_violin(dt, y, hues=[], out_prefix='violin'):
    
    hues = hues if len(hues) else [c for c in dt.columns if c not in ['sample', y]]
    
    for h in hues:
        
        plt.figure(figsize=(2 * len(set(dt[h])), 5))
        
        violin = sns.violinplot(dt, x=h, y=y, hue=h)
        sns.boxplot(dt, x=h, y=y, dodge=False,
                    showcaps=True, boxprops={'facecolor':'None'}, showfliers=False, linewidth=2, whiskerprops={'linewidth':1},
                    ax=violin, legend=False)
        
        #sns.move_legend(violin, loc='upper left', bbox_to_anchor=(1, 1), title=h)
        
        plt.xlabel(None)
        plt.xticks(fontweight='bold')
        
        plt.ylabel(y, fontweight='bold')
        plt.yticks(fontweight='bold')
        
        plt.tight_layout()

        plt.savefig(f'{out_prefix}_{h}_violin.png', dpi=300)
        plt.close()

### ---------------------------------------- ###

def plot_umap(dt, hues=[]):
    
    hues = hues if len(hues) else [c for c in dt.columns if c not in ['sample', 'UMAP_1', 'UMAP_2']]
    
    for h in hues:
        
        legend_cols = int(len(set(dt[h])) / 16) + 1 if h not in ['log10_richness', 'shannon_index'] else 1
        palette = 'tab20' if h not in ['log10_richness', 'shannon_index'] else 'viridis'
        
        plt.figure(figsize=(6, 5))
        
        scatter = sns.scatterplot(dt, x='UMAP_1', y='UMAP_2', hue=h, palette=palette)
        
        sns.move_legend(scatter, loc='upper left', bbox_to_anchor=(1, 1), title=h, ncol=legend_cols)
        
        plt.xlabel('UMAP 1', fontweight='bold')
        plt.xticks(fontweight='bold')
        
        plt.ylabel('UMAP 2', fontweight='bold')
        plt.yticks(fontweight='bold')
        
        plt.tight_layout()

        plt.savefig(f'{h}_umap.png', dpi=300)
        plt.close()

### ---------------------------------------- ###

def plot_pcoa(dt, hues=[]):
    
    hues = hues if len(hues) else [c for c in dt.columns if c not in ['sample', 'PCoA_1', 'PCoA_2', 'UMAP_1', 'UMAP_2']]
    
    for h in hues:
        
        legend_cols = int(len(set(dt[h])) / 16) + 1 if h not in ['log10_richness', 'shannon_index'] else 1
        palette = 'tab20' if h not in ['log10_richness', 'shannon_index'] else 'viridis'
        
        plt.figure(figsize=(6, 5))
        
        scatter = sns.scatterplot(dt, x='PCoA_1', y='PCoA_2', hue=h, palette=palette)
        
        sns.move_legend(scatter, loc='upper left', bbox_to_anchor=(1, 1), title=h, ncol=legend_cols)
        
        plt.xlabel('PCoA 1', fontweight='bold')
        plt.xticks(fontweight='bold')
        
        plt.ylabel('PCoA 2', fontweight='bold')
        plt.yticks(fontweight='bold')
        
        plt.tight_layout()

        plt.savefig(f'{h}_pcoa.png', dpi=300)
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
from scipy.stats import pearsonr
from skbio.stats.ordination import pcoa
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph
from sys import argv

### Parse args

k2_reports, reports_ids, tax_lvl, contaminants, sample_info_path, pca_model_path, umap_model_path, rescale_toggle, use_raw = parse_args()

### Parse reports

print(f'Parsing {len(k2_reports)} reports')

reports_data_raw, reports_data_percent, roots, unclassifieds = [], [], {}, {}
for rid,path in zip(reports_ids, k2_reports):
    
    dt_raw, dt_percent, r, u = parse_k2_report(rid, path, tax_lvl, contaminants)
    
    reports_data_raw.append(dt_raw)
    reports_data_percent.append(dt_percent)
    roots[rid] = r
    unclassifieds[rid] = u

reports_data_raw = pd.concat(reports_data_raw, axis=1)
reports_data_raw.fillna(0, inplace=True)
reports_data_raw.index = reports_data_raw.index.set_names(['clade'])

reports_data_percent = pd.concat(reports_data_percent, axis=1)
reports_data_percent.fillna(0, inplace=True)
reports_data_percent.index = reports_data_percent.index.set_names(['clade'])

if rescale_toggle:

    # Raw
    counts_tot = reports_data_raw.sum(axis=0).values
    norm_factor = np.median(counts_tot) / counts_tot
    reports_data_raw = reports_data_raw.multiply(norm_factor)

    # Percent
    norm_factor = 100 / reports_data_percent.sum(axis=0).values
    reports_data_percent = reports_data_percent.multiply(norm_factor)

reports_data_raw = reports_data_raw.reset_index(drop=False)
reports_data_percent = reports_data_percent.reset_index(drop=False)

reports_data_raw.to_csv('kraken2_summary_raw.tsv', sep='\t', header=True, index=False)
reports_data_percent.to_csv('kraken2_summary_percent.tsv', sep='\t', header=True, index=False)

### Select data for downstream analyses

if use_raw:

    reports_data = reports_data_raw.copy()

else:

    reports_data = reports_data_percent.copy()

### Remove undetected

reports_data = reports_data.loc[reports_data.iloc[:, 1:].sum(axis=1) > 0,].reset_index(drop=True)

### Compute richness and Shannon's index

print('Computing richness and diversity')

stats = {'log10_richness' : [],
         'shannon_index' : []}

for sample in reports_ids:
    
    dt_sub = reports_data.loc[reports_data[sample] > 0, sample].copy()
    
    # Richness
    richness = np.log10(dt_sub.shape[0])
    stats['log10_richness'].append(richness)
    
    # Diversity
    prop = dt_sub.values / 100
    ln_prop = np.log(prop)
    index = - sum(prop * ln_prop)
    stats['shannon_index'].append(index)

stats = pd.DataFrame(stats, index=reports_ids)
stats.index = stats.index.set_names(['sample'])
stats = stats.reset_index(drop=False)

stats.to_csv('sample_stats.tsv', sep='\t', header=True, index=False)

# Plot relationship between Richness and Diversity
r, p = pearsonr(stats.log10_richness.values, stats.shannon_index.values, alternative='two-sided')
m, c = np.polyfit(stats.log10_richness.values, stats.shannon_index.values, 1)

plt.figure(figsize=(5, 5.25))
ax = sns.scatterplot(stats, x='log10_richness', y='shannon_index',
                     size=5, marker='o', color='deepskyblue', edgecolor='black', linewidth=1,
                     legend=False)
fit_line_x = np.linspace(ax.get_xlim()[0],ax.get_xlim()[1],100)
plt.plot(fit_line_x, m * fit_line_x + c, '--', color='black', linewidth=1, alpha=0.75)
plt.xlabel('log$_{10}$(Richness)', fontweight='bold')
plt.ylabel("Shannon's index", fontweight='bold')
plt.title(f'r = {r:.3f}\np = {p:.3e}', loc='left', fontweight='bold')
plt.tight_layout()
plt.savefig('richness_vs_diversity.png', dpi=300)
plt.close()

### PCA + Umap

print('Dimensionality reduction')

pca_data, umap_data = reduce_dimensions(reports_data, pca_model_path, umap_model_path, 50, 30)

pca_data = pca_data.reset_index(drop=False)
umap_data = umap_data.reset_index(drop=False)

pca_data.columns = ['sample'] + pca_data.columns[1:].to_list()
umap_data.columns = ['sample'] + umap_data.columns[1:].to_list()

pca_data.to_csv('pca_coordinates.tsv', sep='\t', header=True, index=False)
umap_data.to_csv('umap_coordinates.tsv', sep='\t', header=True, index=False)

### Pairwise Bray-Curtis dissimilarity

print('Computing pairwise Bray-Curtis dissimilarity')

samples = [c for c in reports_data.columns.values if c != 'clade']

dissimilarity_matrix = np.zeros((len(samples), len(samples)))

for s1 in samples:
    
    for s2 in samples[samples.index(s1) + 1:]:
        
        dissimilarity = 1 - (reports_data[[s1, s2]].min(axis=1).sum() / 100)
        
        dissimilarity_matrix[samples.index(s1), samples.index(s2)] = dissimilarity
        dissimilarity_matrix[samples.index(s2), samples.index(s1)] = dissimilarity

dissimilarity_matrix = pd.DataFrame(dissimilarity_matrix, index=samples, columns=samples)

dissimilarity_matrix.to_csv('dissimilarity_matrix.tsv', sep='\t', index=True, header=True)

### PCoA

print('Generating PCoA embeddings')

pcoa_data = pcoa(dissimilarity_matrix, method='eigh', number_of_dimensions=2, inplace=False).samples

pcoa_data.index = samples
pcoa_data = pcoa_data.reset_index(drop=False)
pcoa_data.columns = ['sample', 'PCoA_1', 'PCoA_2']

pcoa_data.to_csv('pcoa_coordinates.tsv', sep='\t', index=False, header=True)

### Clustering

print('Clustering')

# PCA-based clustering
clusters_pca = cluster_samples(pca_data, 10)
clusters_pca.columns = ['sample', 'pca_clusters']

# Dissimilarity matrix-based clustering
dm = dissimilarity_matrix.reset_index(drop=False)
dm.columns = ['sample'] + dm.columns[1:].to_list()
clusters_dm = cluster_samples(dm, 10)
clusters_dm.columns = ['sample', 'pcoa_clusters']

### Merge datasets

datasets = [pca_data, umap_data, pcoa_data, stats, clusters_pca, clusters_dm]

if sample_info_path != '':
    
    # Load
    sample_info = pd.read_csv(sample_info_path, sep='\t')
    
    datasets.append(sample_info)

try:

    del all_data

except:
    
    pass

for ds in datasets:

    try:

        all_data = pd.merge(all_data, ds, on='sample', how='inner')

    except:
        
        all_data = ds.copy()

### Stats visualizations

print('Generating plots')

hues = [c for c in all_data.columns if c not in ['sample', 'log10_richness', 'shannon_index'] and not c.startswith('PC') and not c.startswith('UMAP')]

# Violin plot for richness
plot_violin(all_data, 'log10_richness', hues=hues, out_prefix='log10_richness')

# Violin plot for diversity
plot_violin(all_data, 'shannon_index', hues=hues, out_prefix='shannon_index')

### Plot UMAP and PCoA data

hues = [c for c in all_data.columns if c != 'sample' and not c.startswith('PC') and not c.startswith('UMAP')]

# UMAP plot
plot_umap(all_data, hues=hues)

# PCoA plot
plot_pcoa(all_data, hues=hues)

### Export updated sample info

all_data.to_csv('updated_sample_info.tsv', sep='\t', index=False, header=True)
