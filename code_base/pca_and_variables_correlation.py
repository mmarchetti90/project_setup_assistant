#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import numpy as np
import pandas as pd
import pickle as pk

from matplotlib import pyplot as plt
from scipy.stats import pearsonr, pointbiserialr, spearmanr
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_classif
from sys import argv

### FUNCTIONS ------------------------------ ###

def parse_args():
    
    # Path to data
    # Must be a comma-delimited or tab-delimited table of shape N * (M + 1)
    # with N = variables, M = samples
    # First column is sample names
    data_path = argv[argv.index('--data') + 1]
    
    # Path to table with variables whose correlation to PCA has to be tested
    # Must be a comma-delimited or tab-delimited table of shape M * (O + 1)
    # with M = samples (same as in data table), O = variables
    # First column is sample names
    # N.B. variable names must have the suffix '_binary', '_categorical' or '_continuous'
    effects_path = argv[argv.index('--effects') + 1]
    
    return data_path, effects_path

### ---------------------------------------- ###

def normalize_data(dat):
    
    # Correct for library size
    dat_sum = dat.sum(axis=0)
    norm_factors = np.median(dat.sum(axis=0)) / dat_sum
    dat = dat * norm_factors
    
    # Log transform
    dat = np.log10(dat + 1)
    
    # Pareto-scaling (mean centering and dividing by sqrt(sd))
    dat = (dat - dat.mean(axis=0)) / np.sqrt(np.std(dat, axis=0))
    
    return dat

### ---------------------------------------- ###

def reduce_dimensions(dat):
    
    # Fit PCA
    
    pca_model = PCA()
    pca_data = pca_model.fit_transform(dat.T)
    
    pca_data = pd.DataFrame(pca_data, index=dat.columns, columns=[f'PC{i+1}' for i in range(pca_data.shape[1])])
    
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
    
    # Export data
    
    pca_data.to_csv('pca_transform.tsv', sep='\t', index=True, header=True)
        
    pk.dump(pca_model, open('pca.pkl', 'wb'))
    
    return pca_model, pca_data, explained_variance, optimal_components

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

def correlate_pc_to_var(pca_dim, var_data):
    
    # Init log
    pc_to_test = pca_dim.columns.values.tolist()
    log, log_header = [], ['var_name', 'var_kind'] + pc_to_test
    
    # Run correlation
    for var in var_data.columns:
        
        # Merge PCA and variable data
        test_data = pd.merge(pca_dim, var_data[var].dropna(), how='inner', left_index=True, right_index=True)
        
        # Test correlation between PC and variable
        if var.endswith('_binary'):
            # Point-Biserial correlation
            binary_vals = np.zeros(test_data.shape[0])
            binary_levels = list(set(test_data[var].dropna().values))
            binary_vals[(test_data[var] == binary_levels[0]).values] = 1            
            pvals = [pointbiserialr(test_data[pc], binary_vals).pvalue for pc in pc_to_test]
            var_name, var_kind = var.replace('_binary', ''), 'binary'
        elif var.endswith('_categorical'):
            # Anova
            pvals = f_classif(test_data[pc_to_test], test_data[var])[1].tolist()
            var_name, var_kind = var.replace('_categorical', ''), 'categorical'
        elif var.endswith('_continuous'):
            # Spearman correlation
            pvals = [spearmanr(test_data[pc], test_data[var]).pvalue for pc in pc_to_test]
            var_name, var_kind = var.replace('_continuous', ''), 'continuous'
        else:
            continue
        log.append([var_name, var_kind] + pvals)
    log = pd.DataFrame(log, columns=log_header)
    
    return log

### MAIN ----------------------------------- ###

if __name__ == '__main__':
    
    # Parse args
    data_path, effects_path = parse_args()
    
    # Load data
    sep = '\t' if data_path.replace('.gz', '').endswith('.tsv') else ','
    data = pd.read_csv(data_path, sep=sep, index_col=0)
    data = data.dropna()
    data.index = data.index.astype(str)
    sep = '\t' if effects_path.replace('.gz', '').endswith('.tsv') else ','
    effects = pd.read_csv(effects_path, sep=sep, index_col=0)
    effects.index = effects.index.astype(str)
    
    # Normalize data
    norm_data = normalize_data(data)
    
    # PCA transform
    pca_model, pca_data, explained_variance, optimal_components = reduce_dimensions(norm_data)
    pc_to_test = [f'PC{i+1}' for i in range(optimal_components + 1)]
    
    # Correlate PCA dimensions to variables
    var_correl = correlate_pc_to_var(pca_data[pc_to_test], effects)
    
    # Export reports
    var_correl.to_csv('pca_correlation.tsv', sep='\t', index=False)        
