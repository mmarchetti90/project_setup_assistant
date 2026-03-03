#!/usr/bin/env python3

"""
This script runs a Differential Expression Analysis using U Test
N.B. Using mannwhitneyu instead of ranksums since the latter does not handle ties
N.B. Only use if sample size > 8 for each sample type

INPUTS
------

counts_file:
    Table N * (M + 1), where N is the number of genes, and m is the number of columns
    The first column should be named GeneID and report unique identifiers for genes
    
comparison_file:
    Tab-separated file with the following structure:
        ComparisonName      <output_files_prefix>
        Reference           <reference_condition>
        sample              condition
        <sample_1_id>       <condition>
        ...                 ...
        <sample_M_id>       <condition>
"""

### ---------------------------------------- ###

def parse_args():
    
    # Load counts matrix
    counts_file_path = argv[argv.index("--counts_file") + 1]
    counts = pd.read_csv(counts_file_path, sep='\t', header=0)
    
    # Read comparison file
    comparison_file_path = argv[argv.index("--comparison_file") + 1]
    comparison_info = parse_comparison_file(comparison_file_path)

    # Remove files from comparison_info not in counts table
    comparison_info[1]['samples'] = [sample for sample in comparison_info[1]['samples'] if sample in counts.columns.to_list()]
    comparison_info[1]['reference'] = [sample for sample in comparison_info[1]['reference'] if sample in counts.columns.to_list()]
    
    # Remove files from counts table not in comparison_info
    counts = counts.loc[:, ['GeneID'] + comparison_info[1]['samples'] + comparison_info[1]['reference']]
    
    # Min reads per gene across samples
    if '--min_reads' in argv:
        
        min_reads = int(argv[argv.index("--min_reads") + 1])
    
    else:
        
        min_reads = 10
    
    # Filter genes
    counts = counts.loc[(counts.iloc[:, 1:] >= min_reads).sum(axis=1) > (counts.shape[1] / 2), ]
    
    # Significance threshold
    if '--significance_threshold' in argv:
        
        p_thr = float(argv[argv.index("--significance_threshold") + 1])
    
    else:
        
        p_thr = 0.05

    if '--norm_mode' in argv:

        norm_mode = argv[argv.index('--norm_mode') + 1]

    else:

        norm_mode = 'PFlog1pPF'
        
    return counts, comparison_info, p_thr, norm_mode

### ---------------------------------------- ###

def parse_comparison_file(comparison_file_path):
    
    raw_info = open(comparison_file_path).read().split('\n')
    
    comparison_name = raw_info[0].split('\t')[1]
    
    reference = raw_info[1].split('\t')[1]
    
    samples = {'samples' : [], 'reference' : []}
    
    for line in raw_info[3:]:
        
        sample_id, sample_type = line.split('\t')[:2]
        
        if sample_type == reference:
            
            samples['reference'].append(sample_id)
        
        else:
            
            samples['samples'].append(sample_id)
    
    return [comparison_name, samples]

### ---------------------------------------- ###

def normalize_counts(raw_counts, mode='PFlog1pPF'):
    
    if mode not in ['DESeq2', 'PFlog1pPF']:
        
        print('ERROR: unknown normalization method.\nAccepted values are "DESeq2" and "PFlog1pPF"')
        normalized_counts = pd.DataFrame()
    
    elif mode == 'DESeq2':
        
        normalized_counts = raw_counts.copy().replace(0, 1, inplace=False)
        
        normalized_counts = normalized_counts.astype({col : (str if n == 0 else float) for n,col in enumerate(normalized_counts.columns.values)})
        
        # Deseq2 normalization
        pseudo_ref_sample = gmean(normalized_counts.iloc[:, 1:].replace(0, 1, inplace=False), axis=1)
        norm_factors = np.median(normalized_counts.iloc[:, 1:].div(pseudo_ref_sample, axis=0), axis=0)
        normalized_counts.iloc[:, 1:] = normalized_counts.iloc[:, 1:].div(norm_factors, axis=1)
    
    else:
        
        normalized_counts = raw_counts.copy()
        
        normalized_counts = normalized_counts.astype({col : (str if n == 0 else float) for n,col in enumerate(normalized_counts.columns.values)})
        
        # Proportional fitting
        library_sizes = normalized_counts.iloc[:, 1:].sum(axis=0).values
        median_size = np.median(library_sizes)
        normalized_counts.iloc[:, 1:] = normalized_counts.iloc[:, 1:].div(library_sizes / median_size, axis=1)
        
        # log1p
        normalized_counts.iloc[:, 1:] = np.log1p(normalized_counts.iloc[:, 1:])
        
        # Proportional fitting
        library_sizes = normalized_counts.iloc[:, 1:].sum(axis=0).values
        median_size = np.median(library_sizes)
        normalized_counts.iloc[:, 1:] = normalized_counts.iloc[:, 1:].div(library_sizes / median_size, axis=1)
    
    return normalized_counts

### ---------------------------------------- ###

def differential_expression(analysis_name, cnts, sample_ids, ctrl_ids, p_threshold=0.05, norm_mode='PFlog1pPF'):

    # Data normalization
    cnts = normalize_counts(cnts, mode=norm_mode)

    # Differential expression analysis
    sample_mean = cnts.loc[:, sample_ids].mean(axis=1)
    
    ctrl_mean = cnts.loc[:, ctrl_ids].mean(axis=1)

    sample_median = cnts.loc[:, sample_ids].median(axis=1)
    
    ctrl_median = cnts.loc[:, ctrl_ids].median(axis=1)
    
    if norm_mode == 'PFlog1pPF':

        log2fc = sample_median - ctrl_median

    else:

        log2fc = np.log2(sample_median / ctrl_median)
    
    pval = mannwhitneyu(cnts.loc[:, sample_ids], cnts.loc[:, ctrl_ids], axis=1).pvalue

    fdr = fdrcorrection(pval, alpha=0.05, is_sorted=False)[1]
    
    results = pd.DataFrame({'sample_mean' : sample_mean,
                            'ctrl_mean' : ctrl_mean,
                            'sample_median' : sample_median,
                            'ctrl_median' : ctrl_median,
                            'log2FoldChange' : log2fc,
                            'pvalue' : pval,
                            'padj' : fdr})

    results.index = cnts.GeneID.values
    
    results.to_csv(f'DEA_{analysis_name}.tsv', sep='\t', header=True, index=True)
    
    # Sample correlation plot
    #correlation_analysis(f'{analysis_name}_CorrelationPlot.png', cnts.loc[:, sample_ids + ctrl_ids])
    
    # PCA
    pca_analysis(f'{analysis_name}_PCA.png', cnts.loc[:, sample_ids + ctrl_ids].T, len(sample_ids), len(ctrl_ids))
    
    # MA plot
    ma_plot(f'{analysis_name}_MA-Plot.png', 0.5 * (np.log2(results.sample_mean.values + 1) + np.log2(results.ctrl_mean.values + 1)), results.log2FoldChange.values)
    
    # Variance vs mean plot
    #variance_vs_mean_plot(f'{analysis_name}_VarianceVsMean.png', cnts.iloc[:, 1:].mean(axis=0), cnts.iloc[:, 1:].var(axis=0))
    
    # Volcano plot
    volcano_plot(f'{analysis_name}_VolcanoPlot.png', results, p_threshold)
    
### ---------------------------------------- ###

def correlation_analysis(out_name, data):
    
    sns.set_theme(style="white")

    # Compute the correlation matrix
    corr = data.corr(method='pearson')
    
    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(0.5 * corr.shape[0], 0.5 * corr.shape[0]))
    
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-1, vmax=1, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})
    
    plt.xticks(rotation=45)
    
    plt.savefig(out_name, bbox_inches='tight', dpi=300)
    plt.close()

### ---------------------------------------- ###

def pca_analysis(out_name, data, n_samples, n_ctrls):
    
    pca_model = PCA().fit(data)
    pca_data = pca_model.transform(data)
    
    explained_variance_1, explained_variance_2 = pca_model.explained_variance_[:2] * (100 / pca_model.explained_variance_.sum())
    
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=pca_data[:, 0], y=pca_data[:, 1], hue=[out_name.split('_')[0] for _ in range(n_samples)] + [out_name.split('_')[2] for _ in range(n_ctrls)])
    plt.xlabel(f'PC1\n{explained_variance_1}%')
    plt.ylabel(f'PC2\n{explained_variance_2}%')
    plt.savefig(out_name, bbox_inches='tight', dpi=300)
    plt.close()

### ---------------------------------------- ###

def ma_plot(out_name, log_mean, log_fc):
    
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=log_mean, y=log_fc)
    plt.xlabel('Log2FC')
    plt.ylabel('Mean of normalized counts')
    plt.savefig(out_name, bbox_inches='tight', dpi=300)
    plt.close()

### ---------------------------------------- ###

def variance_vs_mean_plot(out_name, gene_mean, gene_var):
    
    plt.figure(figsize=(10, 10))
    sns.scatterplot(x=gene_mean, y=gene_var)
    plt.xlabel('Gene mean')
    plt.ylabel('Gene variance')
    plt.savefig(out_name, bbox_inches='tight', dpi=300)
    plt.close()

### ---------------------------------------- ###

def volcano_plot(out_name, data, p_threshold):
    
    # Get number of upregulated and downregulated genes
    n_up = data.loc[(data.log2FoldChange > 0) & (data.padj < p_threshold),].shape[0]
    n_down = data.loc[(data.log2FoldChange < 0) & (data.padj < p_threshold),].shape[0]
    n_ns = data.shape[0] - n_up - n_down
    
    # Add useful columns
    data['transformed_padj'] = - np.log10(data.padj.values)
    data['DEA'] = [f'ns ({n_ns})' for _ in range(data.shape[0])]
    data.loc[(data.log2FoldChange > 0) & (data.padj < p_threshold), 'DEA'] = f'upreg ({n_up})'
    data.loc[(data.log2FoldChange < 0) & (data.padj < p_threshold), 'DEA'] = f'downreg ({n_down})'
    
    plt.figure(figsize=(10, 10))
    sns.scatterplot(data=data, x='log2FoldChange', y='transformed_padj', hue='DEA', palette=['red', 'lightgray', 'green'], hue_order=[f'downreg ({n_down})', f'ns ({n_ns})', f'upreg ({n_up})'], linewidth=0.5, edgecolor='black')
    plt.xlabel('Log2FC')
    plt.ylabel('- Log10 pvalue')
    plt.savefig(out_name, bbox_inches='tight', dpi=300)
    plt.close()

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.stats import gmean, mannwhitneyu
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import fdrcorrection
from sys import argv

# Import data
counts, comparison_info, p_thr, norm_mode = parse_args()

# Differential expression analysis
differential_expression(comparison_info[0], counts, comparison_info[1]['samples'], comparison_info[1]['reference'], p_threshold=p_thr, norm_mode=norm_mode)
