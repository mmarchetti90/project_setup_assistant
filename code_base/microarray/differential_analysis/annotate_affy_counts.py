#!/usr/bin/env python3

"""
This script annotates affy counts
"""

### ---------------------------------------- ###

def parse_args():
    
    count_files = argv[argv.index('--count_files') + 1].split(';')

    gtf_file = argv[argv.index('--gtf') + 1]
    
    if '--chromosomes' in argv:
        
        chromosomes_of_interest = argv[argv.index('--chromosomes') + 1].split(',')
    
    else:
        
        chromosomes_of_interest = []
        
    sample_info_file = argv[argv.index('--info') + 1]
    
    return count_files, gtf_file, chromosomes_of_interest, sample_info_file

### ---------------------------------------- ###

def parse_gtf(gtf_file, chromosomes_of_interest=[]):
    
    raw_data = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, dtype=str)
    raw_data.loc[:, 0] = raw_data.loc[:, 0].astype(str)
    raw_data = raw_data.loc[raw_data.iloc[:, 2] == 'gene',]
    
    genes_info = {'gene_id' : [],
                  'gene_symbol' : [],
                  'chromosome' : [],
                  'biotype' : []}
    
    for n,entry in raw_data.iterrows():
        
        chrom, _, _, _, _, _, _, _, info = entry.values
        
        if len(chromosomes_of_interest) > 0 and chrom not in chromosomes_of_interest:
            
            continue
        
        try:
            
            gene_id = [txt[len('gene_id "') : -1].replace(';', '').replace('"', '') for txt in info.split('; ') if 'gene_id' in txt][0]
            
        except:
            
            gene_id = ''
    
        
        try:
            
            gene_symbol = [txt[len('gene_name "') : -1].replace(';', '').replace('"', '') for txt in info.split('; ') if 'gene_name' in txt][0]
            
        except:
            
            gene_symbol = ''
            
        
        try:
            
            biotype = [txt[len('gene_biotype "') : -1].replace(';', '').replace('"', '') for txt in info.split('; ') if 'gene_biotype' in txt][0]
            
        except:
            
            biotype = ''
        
        genes_info['gene_id'].append(gene_id)
        genes_info['gene_symbol'].append(gene_symbol)
        genes_info['chromosome'].append(chrom)
        genes_info['biotype'].append(biotype)
    
    genes_info = pd.DataFrame(genes_info)
    
    return genes_info

### ------------------MAIN------------------ ###

import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from sys import argv

### Parse CLI args

count_files, gtf_file, chromosomes_of_interest, sample_info_file = parse_args()

### Load GTF and info files

gtf = parse_gtf(gtf_file, chromosomes_of_interest)
sample_info = pd.read_excel(sample_info_file, comment='#')

### Load and merge count files (columns with same names are averaged)

# Merge data
try:
    
    del all_data
    
except:
    
    pass


for file in count_files:
    
    new_data = pd.read_csv(file, sep='\t', index_col=0, header=0)
    
    try:
        
        all_data = all_data.merge(new_data, left_index=True, right_index=True, how='inner')
    
    except:
        
        all_data = new_data.copy()

# Correct column names, then group and average
all_data.columns = [col.replace('_x', '').replace('_y', '') for col in all_data.columns]
all_data = all_data.groupby(by=all_data.columns, axis=1).mean()

### Add gene_names, biotypes, and chromosome

# Subset gtf to remove genes not in all_data, then do the opposite with all_data, reaordering the latter based on gtf genes order
gtf = gtf.loc[gtf.gene_symbol.isin(all_data.index),]
all_data = all_data.loc[gtf.gene_symbol.values,]

# Add columns for gene_symbol, gene_biotype, and chromosome
og_columns = all_data.columns.to_list()
all_data = all_data.assign(gene_id = gtf.gene_id.values,
                           gene_symbol = gtf.gene_symbol.values,
                           gene_biotype = gtf.biotype.values,
                           chromosome = gtf.chromosome.values)
all_data = all_data.loc[:, ['gene_id', 'gene_symbol', 'gene_biotype', 'chromosome'] + og_columns]

# sum(all_data.index != all_data.gene_id) # debug

### Add sample info

sample_type = [sample_info.loc[sample_info.File == col.replace('.gz', ''), 'Type'].values[0] if col.replace('.gz', '') in sample_info.File.values else '' for col in all_data.columns]
sample_type = pd.DataFrame([sample_type], columns=all_data.columns)

trisomy_status = [sample_info.loc[sample_info.File == col.replace('.gz', ''), 'tri(8)'].values[0] if col.replace('.gz', '') in sample_info.File.values else '' for col in all_data.columns]
trisomy_status = ['tri(8)' if ts == 1 else '' if ts == '' else 'Normal' for ts in trisomy_status]
trisomy_status = pd.DataFrame([trisomy_status], columns=all_data.columns)

all_data = pd.concat([sample_type, trisomy_status, all_data], ignore_index=True)

### Sort columns by sample type and trisomy status

columns_order = (['gene_id', 'gene_symbol', 'gene_biotype', 'chromosome'] +
                 [col for col,st,ts in zip(all_data.columns, sample_type.values[0], trisomy_status.values[0]) if st == 'Normal'] +
                 [col for col,st,ts in zip(all_data.columns, sample_type.values[0], trisomy_status.values[0]) if st == 'MDS' and ts == 'Normal'] +
                 [col for col,st,ts in zip(all_data.columns, sample_type.values[0], trisomy_status.values[0]) if st == 'MDS' and ts != 'Normal'])

all_data = all_data.loc[:, columns_order]

### Save to file

all_data.to_csv('annotated_count_data.tsv', sep='\t', index=False, header=True)
