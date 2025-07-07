#!/usr/bin/env python3

"""
Scoring gene signature of interest in various subpopulations
Note that this is intended for large sample sizes

INPUT
-----

gene_counts
    Tab-separated matrix of raw gene counts with shape N * (M + 1) or N * (M + 2)
    One of two columns (or both) must be present as leftmost column(s): GeneSymbol and GeneID

gene_list
    Plain text file with one gene per row

gene_list_name
    Output prefix

sample_info
    Optional tab-separated matrix with two columns with shape N' * (M + 1), with N' being a subset
    of N samples and M the columns with samples groupings
    The first column reports samples' names (matching gene_signature)
"""

### ---------------------------------------- ###

def parse_args():

    print(argv)
    
    # Load gene counts
    gene_counts_path = argv[argv.index('--gene_counts') + 1]
    gene_counts = pd.read_csv(gene_counts_path, sep='\t')
    
    # Modify TCGA-COAD gene ids
    gene_counts['GeneID'] = [val.split('.')[0] for val in gene_counts['GeneID'].values]
    
    # Add GeneID or GeneSymbol columns if one is missing (makes a copy of the available one)
    if 'GeneID' in gene_counts.columns and 'GeneSymbol' in gene_counts.columns:
        
        pass
    
    elif 'GeneID' not in gene_counts.columns:
        
        gene_counts['GeneID'] = gene_counts['GeneSymbol'].values
    
    else:
        
        gene_counts['GeneSymbol'] = gene_counts['GeneID'].values
        
    gene_counts = gene_counts.loc[:, ['GeneID', 'GeneSymbol'] + [col for col in gene_counts.columns if col not in ['GeneID', 'GeneSymbol']]]
    
    # Load 
    gene_list_path = argv[argv.index('--gene_list') + 1]
    gene_list = np.array([gene for gene in open(gene_list_path, 'r').read().split('\n') if len(gene)])
    
    # Remove genes not found in the counts matrix
    if 'GeneID' in gene_counts.columns:
        
        filter_1 = np.isin(gene_list, gene_counts['GeneID'].values, assume_unique=True)
    
    else:
        
        filter_1 = np.zeros(len(gene_list))
    
    if 'GeneSymbol' in gene_counts.columns:
        
        filter_2 = np.isin(gene_list, gene_counts['GeneSymbol'].values, assume_unique=True)
    
    else:
        
        filter_2 = np.zeros(len(gene_list))
    
    print(f'Found {(filter_1 | filter_2).sum()} / {len(gene_list)} genes from list')
    
    gene_list = gene_list[filter_1 | filter_2]
    
    gene_list = np.unique(gene_list)
    
    # Gene list name (used for outputs)
    gene_list_name = argv[argv.index('--gene_list_name') + 1].replace(' ', '_')
    
    # Sample info
    if '--sample_info' in argv:

        sample_info_path = argv[argv.index('--sample_info') + 1]
        sample_info = pd.read_csv(sample_info_path, sep='\t', header=0, index_col=0)
    
    else:
        
        sample_info = pd.DataFrame({})
    
    return gene_counts, gene_list, gene_list_name, sample_info

### ---------------------------------------- ###

def normalize_counts(raw_counts):
    
    # Copying raw_counts matrix
    normalized_counts = raw_counts.copy().astype({col : (str if col in ['GeneID', 'GeneSymbol'] else float) for col in raw_counts.columns})
    
    # Proportional fitting
    counts_per_cell = np.array(normalized_counts.iloc[:, 2:].sum(axis = 0))
    median_counts = np.median(counts_per_cell)
    normalized_counts.iloc[:, 2:] = normalized_counts.iloc[:, 2:].div(counts_per_cell / median_counts, axis = 1)
    
    # log1p
    normalized_counts.iloc[:, 2:] = np.log1p(normalized_counts.iloc[:, 2:])
    
    # Proportional fitting
    counts_per_cell = np.array(normalized_counts.iloc[:, 2:].sum(axis = 0))
    median_counts = np.median(counts_per_cell)
    normalized_counts.iloc[:, 2:] = normalized_counts.iloc[:, 2:].div(counts_per_cell / median_counts, axis = 1)
    
    return normalized_counts

### ---------------------------------------- ###

def scale_features(normalized_counts):
    
    # Copying normalized_counts matrix
    scaled_features = normalized_counts.copy()
    
    # Scaling features
    features_mean = normalized_counts.iloc[:, 2:].mean(axis = 1).to_numpy()
    features_std = normalized_counts.iloc[:, 2:].std(axis = 1).to_numpy()
    scaled_features.iloc[:, 2:] = scaled_features.iloc[:, 2:].subtract(features_mean, axis="index").div(features_std, axis="index")
    
    return scaled_features

### ---------------------------------------- ###

def score_gene_set(scaled_features, gene_set, gene_pool=[], bins=50, ctrl_genes_num=50):
    
    # Reset index
    scaled_features = scaled_features.reset_index(drop=True)

    # Subsetting gene_set and gene_list for detected genes
    gene_set = [gene for gene in gene_set if gene in scaled_features.GeneID.tolist() or gene in scaled_features.GeneSymbol.tolist()]
    gene_pool = [gene for gene in gene_pool if gene in scaled_features.GeneID.tolist() or gene in scaled_features.GeneSymbol.tolist()]

    if not len(gene_pool):
        
        gene_pool = scaled_features.loc[(~ scaled_features.GeneID.isin(gene_set)) &
                                        (~ scaled_features.GeneSymbol.isin(gene_set)), 'GeneID'].to_list()
    
    # Mean of each gene in the gene_pool across all cells
    gene_means = scaled_features.loc[(scaled_features.GeneID.isin(gene_pool + gene_set)) |
                                     (scaled_features.GeneSymbol.isin(gene_pool + gene_set)),].iloc[:, 2:].mean(axis=1)
    
    # Rank genes based on binned expression level, then for each bin of the genes in the gene_set, pick ctrl_genes random genes for genes with matched binned expression
    bin_size = len(gene_means) // bins
    expression_order = np.argsort(gene_means)
    
    set_indexes = expression_order.loc[(scaled_features.GeneID.isin(gene_set)) |
                                       (scaled_features.GeneSymbol.isin(gene_set))].values.tolist()
    
    ctrl_set = []
    for index in set_indexes:
        
        random_pool = expression_order.index[(expression_order >= (index - bin_size // 2)) &
                                             (expression_order <= (index + bin_size // 2)) &
                                             ~(expression_order.isin(set_indexes))].tolist()
        np.random.shuffle(random_pool)
        ctrl_set.extend(random_pool[:ctrl_genes_num])
    ctrl_set = list(set(ctrl_set)) # Removing duplicates
    
    ctrl_set = scaled_features.loc[ctrl_set, 'GeneID']
    
    # Computing the mean of gene_set and ctrl_set genes for each cell
    set_means = scaled_features.loc[(scaled_features.GeneID.isin(gene_set)) |
                                    (scaled_features.GeneSymbol.isin(gene_set)),].iloc[:, 2:].mean(axis=0)
    ctrl_means = scaled_features.loc[(scaled_features.GeneID.isin(ctrl_set)) |
                                     (scaled_features.GeneSymbol.isin(ctrl_set)),].iloc[:, 2:].mean(axis=0)
    
    set_score = set_means - ctrl_means
    
    return set_score

### ---------------------------------------- ###

def get_stats(c1dt, c2dt):
    
    # Stats
    c1_mean, c1_median, c1_std = np.mean(c1dt), np.median(c1dt), np.std(c1dt)
    c2_mean, c2_median, c2_std = np.mean(c2dt), np.median(c2dt), np.std(c2dt)
    means_delta = c1_mean - c2_mean
    pval = mannwhitneyu(c1dt, c2dt, alternative='two-sided').pvalue
    
    return c1_mean, c1_median, c1_std, c2_mean, c2_median, c2_std, means_delta, pval

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from sys import argv

### Parse args

gene_counts, gene_list, gene_list_name, sample_info = parse_args()

### Pre-process counts

# Remove genes not detected in any sample
gene_counts = gene_counts.loc[gene_counts.iloc[:, 2:].sum(axis=1) != 0, ]

# Normalization (PFlog1pPF)
norm_counts = normalize_counts(gene_counts)

# Scale data
scaled_counts = scale_features(norm_counts)

### Calculate gene signatures

gene_signature = score_gene_set(scaled_counts, gene_list, bins=50, ctrl_genes_num=(20 * len(gene_list.shape)))

gene_signature.name = gene_list_name

### Save to file

gene_signature.to_csv(f'{gene_list_name}.tsv', sep='\t', index=True, header=True)

### Filter and merge data

data = pd.merge(sample_info, gene_signature, left_index=True, right_index=True, how='inner')

### Stats

stats_keys = ['condition_1', 'condition_2',
              'c1_mean', 'c1_median', 'c1_std',
              'c2_mean', 'c2_median', 'c2_std',
              'means_delta', 'pval', 'padj']

for col in data.columns[:-1]:
    
    out_prefix = f'{gene_list_name}_{col}'
    
    # Init stats dict
    gene_signature_stats = {k : [] for k in stats_keys}
    p_subset = []
    
    # Stats
    unique_conditions = list(dict.fromkeys(data[col].to_list()).keys()) # This preserves the order in sample_info
    for c1 in unique_conditions[:-1]:
        
        for c2 in unique_conditions[unique_conditions.index(c1) + 1:]:
            
            # Subset data
            c1_dat = data.loc[data[col].values == c1, gene_list_name].values
            c2_dat = data.loc[data[col].values == c2, gene_list_name].values
            
            # Analyze
            stats_values = get_stats(c1_dat, c2_dat)
            
            # Store stats data
            for k,v in zip(stats_keys, [c1, c2] + list(stats_values)):
                
                gene_signature_stats[k].append(v)
                
                if k == 'pval':
                    
                    p_subset.append(v)
            
    gene_signature_stats['padj'] = fdrcorrection(p_subset, alpha=0.05, is_sorted=False)[1]
    
    gene_signature_stats = pd.DataFrame(gene_signature_stats)

    gene_signature_stats.to_csv(f'{out_prefix}_stats.tsv', sep='\t', index=False)
    
### Plotting

for col in data.columns[:-1]:
    
    out_prefix = f'{gene_list_name}_{col}'
    
    hue_order = list(dict.fromkeys(data[col].to_list()).keys()) # This preserves the order in sample_info
        
    plt.figure(figsize=(15, 5))

    #ax = sns.stripplot(data, x=col, y=signature_name, hue=col, dodge=False, hue_order=hue_order)
    ax = sns.swarmplot(data, x=col, y=gene_list_name, hue=col, dodge=False, hue_order=hue_order)
    sns.boxplot(data, x=col, y=gene_list_name,
                showcaps=False, boxprops={'facecolor':'None'}, showfliers=False, whiskerprops={'linewidth':0}, ax=ax)

    plt.xlabel(None)
    plt.ylabel('Signature', c='black', fontweight='bold', fontsize=15)
    plt.xticks(c='black', fontweight='bold', fontsize=15)

    plt.tight_layout()
    plt.savefig(f'{out_prefix}.png', dpi=300)
    plt.close()
