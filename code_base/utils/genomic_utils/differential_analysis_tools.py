#!/usr/bin/env python3

"""
Collection of utilities for parsing Differential Analysis (DA) results
(e.g. expression, accessibility, etc)
"""

### ---------------------------------------- ###

def load_and_split(path: str, fc_col: str, p_col: str, fc_thr: float=0., p_thr: float=0.05):
    
    """
    Function that loads a DA result file, filters data based on adjusted p value and splits data on
    a fold change threshold
    
    Parameters
    ----------
    path : str
        Path to DA table
    fc_col : str
        Log2 Fold Change column
    padj_col : str
        Adjusted p value column
    fc_thr : float
        Log2 fold change absolute threshold
    p_thr : float
        Adjusted p value threshold
    """
    
    data = pd.read_csv(path, sep='\t')
    
    upreg = data.loc[(data[p_col] < p_thr) & (data[fc_col] > fc_thr),].copy()
    
    downreg = data.loc[(data[p_col] < p_thr) & (data[fc_col] < - fc_thr),].copy()
    
    return upreg, downreg

### ---------------------------------------- ###

def load_and_merge(paths: dict, common_cols: list, fc_col: str, padj_col: str):
    
    """
    Function that loads DA results files and merges them based on a set of common columns
    (e.g. gene name, gene symbol, etc)
    
    Parameters
    ----------
    paths : dict
        Dictionary with unique sample identifiers as keys and paths as values
    common_cols : list
        Common columns to use for merging datasets
    fc_col : str
        Log2 Fold Change column
    padj_col : str
        Adjusted p value column
    """

    try:
        
        del all_data
        
    except:
        
        pass
    
    for sample_name,path in paths.items():
    
        data = pd.read_csv(path, sep='\t')
        
        data = data[common_cols + [fc_col, padj_col]]
        
        data.columns = common_cols + [f'{fc_col}_{sample_name}', f'{padj_col}_{sample_name}']
        
        try:
            
            all_data = pd.merge(all_data, data.copy(), on=common_cols, how='outer')
        
        except:
            
            all_data = data.copy()
    
    all_data[[col for col in all_data.columns if col.startswith(padj_col)]] = all_data[[col for col in all_data.columns if col.startswith(padj_col)]].fillna(1)

    all_data[[col for col in all_data.columns if col.startswith(fc_col)]] = all_data[[col for col in all_data.columns if col.startswith(fc_col)]].fillna(0)
    
    return all_data

### ---------------------------------------- ###

try:
    
    import pandas as pd
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()
