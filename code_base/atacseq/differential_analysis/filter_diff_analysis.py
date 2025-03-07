#!/usr/bin/env python3

"""
This script filters ATAC differential analysis results by:
- removing genes with peaks both more and less open in the same condition (i.e. unclear net effect);
- merging peaks for the same gene that have the same behaviour (i.e. net positive/negative effect, largest effect will be kept);
"""

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

import pandas as pd

from sys import argv

### Load differential analysis

p_thr = 0.01

print(argv)
diff_file = argv[argv.index('--differential_analysis') + 1]
diff_data = pd.read_csv(diff_file, sep='\t', header=0)
diff_data.loc[diff_data.GeneName.isna(), 'GeneName'] = diff_data.loc[diff_data.GeneName.isna(), 'GeneID'] # Replace missing GeneNames with GeneIDs

### Fill NAs in padj column

diff_data.loc[diff_data['padj'].isna(), 'padj'] = 1

### Remove peaks without assigned genes

diff_data = diff_data.loc[~ pd.isna(diff_data.GeneID),]

### Filter genes based on peaks

rows_to_be_removed = []
genes_one_peak, genes_multi_peak, genes_discordant_peak = [], [], []
for gene in set(diff_data.GeneID):
    
    gene_peaks = diff_data.loc[diff_data.GeneID == gene,].copy()
    significant_gene_peaks = gene_peaks.loc[gene_peaks.padj < p_thr,].copy()
    
    if gene_peaks.shape[0] == 1: # Only one peak
    
        if gene_peaks.padj.values[0] < p_thr:
            
            genes_one_peak.append(gene)
    
    elif significant_gene_peaks.shape[0] == 1: # Only one significant peak
        
        best_row = significant_gene_peaks.index[0]
        
        bad_rows = [i for i in gene_peaks.index if i != best_row]
        rows_to_be_removed.extend(bad_rows)
        
        genes_one_peak.append(gene)
    
    elif significant_gene_peaks.shape[0] == 0:
        
        if sum(gene_peaks.log2FoldChange.values > 0) == gene_peaks.shape[0] or sum(gene_peaks.log2FoldChange.values < 0) == gene_peaks.shape[0]:
            
            best_row = gene_peaks.index[gene_peaks.padj == gene_peaks.padj.min()][0]
            
            bad_rows = gene_peaks.index[gene_peaks.padj != gene_peaks.padj.min()].to_list()
            rows_to_be_removed.extend(bad_rows)
        
        else:
            
            bad_rows = gene_peaks.index.to_list()
            rows_to_be_removed.extend(bad_rows)
    
    else:
        
        if sum(significant_gene_peaks.log2FoldChange.values > 0) == significant_gene_peaks.shape[0] or sum(significant_gene_peaks.log2FoldChange.values < 0) == significant_gene_peaks.shape[0]:
            
            best_row = significant_gene_peaks.index[significant_gene_peaks.log2FoldChange.abs() == significant_gene_peaks.log2FoldChange.abs().max()][0]
            
            bad_rows = gene_peaks.index[gene_peaks.log2FoldChange.abs() != gene_peaks.log2FoldChange.abs().max()].to_list()
            rows_to_be_removed.extend(bad_rows)
            
            genes_multi_peak.append(gene)
        
        else:
            
            bad_rows = gene_peaks.index.to_list()
            rows_to_be_removed.extend(bad_rows)
            
            genes_discordant_peak.append(gene)

largest_cat = max([len(genes_one_peak), len(genes_multi_peak), len(genes_discordant_peak)])
genes_categories = pd.DataFrame({"one_peak" : genes_one_peak + ['' for _ in range(largest_cat - len(genes_one_peak))],
                                 "multi_peak_concordant" : genes_multi_peak + ['' for _ in range(largest_cat - len(genes_multi_peak))],
                                 "multi_peak_discordant" : genes_discordant_peak + ['' for _ in range(largest_cat - len(genes_discordant_peak))]})
genes_categories.to_csv('genes_categories.tsv', header=True, index=False, sep='\t')

diff_data.drop(index=rows_to_be_removed, inplace=True)

### Remove genes with same name (small RNAs, etc)

diff_data = diff_data.loc[~ diff_data.GeneName.duplicated(keep=False),]

### Export

diff_data.to_csv("diff_analysis_filtered.tsv", index=False, sep='\t')
