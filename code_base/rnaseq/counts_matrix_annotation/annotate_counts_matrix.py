#!/usr/bin/env python3

"""
This script annotates a raw counts matrix with gene ID, symbol, biotype, and coordinates

Input counts table should be tab-separated with the first column reporting unique gene identifiers
and the following should report the gene counts for individual samples
"""

### ---------------------------------------- ###

def load_gtf(path, desired_biotypes=[], desired_chromosomes=[]):
    
    # Load GTF
    
    gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    gtf_data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Only keep genes

    gtf_data = gtf_data.loc[gtf_data.feature == 'gene', ['seqname', 'start', 'end', 'strand', 'attribute']]
    
    # Get biotype and gene id

    biotypes, gene_ids, gene_symbols = [], [], []
    for _,row in gtf_data.iterrows():
        
        info = row.values[-1]
        
        biotype = re.findall('gene_biotype "\w+";', info)[0]
        biotype = biotype.replace('gene_biotype ', '').replace(';', '').replace('"', '')
        
        biotypes.append(biotype)
        
        gene = re.findall('gene_id "\w+";', info)[0]
        gene = gene.replace('gene_id ', '').replace(';', '').replace('"', '')
        
        gene_ids.append(gene)
        
        if 'gene_name' in info:
            
            gene = info[info.index('gene_name "') + len('gene_name "'):]
            gene = gene[:gene.index('"')]
        
        else:
            
            gene = ''
        
        gene_symbols.append(gene)

    gtf_data['biotype'] = biotypes
    gtf_data['gene_id'] = gene_ids
    gtf_data['gene_symbol'] = gene_symbols
    
    # Filter based on biotype
    
    if len(desired_biotypes):

        gtf_data = gtf_data.loc[gtf_data.biotype.isin(desired_biotypes),]

    # Filter for desired chromosomes

    if len(desired_chromosomes):

        gtf_data = gtf_data.loc[gtf_data.seqname.isin(desired_chromosomes),]
    
    # Remove genes without gene_symbol
    
    #gtf_data = gtf_data.loc[gtf_data['gene_symbol'] != '',]
    
    # Fix dtypes
    
    gtf_data[['start', 'end']] = gtf_data[['start', 'end']].astype(int)
    
    return gtf_data

### ------------------MAIN------------------ ###

import pandas as pd
import re

from sys import argv

### Parse arguments

counts_file = argv[argv.index("--counts_path") + 1]
gtf_file = argv[argv.index("--gtf_path") + 1]

### Import counts file

counts = pd.read_csv(counts_file, sep='\t', header=0)

counts.columns = ['gene_id'] + counts.columns[1:].to_list()

### Parse gtf file to find gene sizes (sum of exon lengths)

gtf = load_gtf(gtf_file)

gtf = gtf[['gene_id', 'gene_symbol', 'biotype', 'seqname', 'start', 'end', 'strand']]

### Merge and save annotated counts to file

merged_data = pd.merge(gtf, counts, how='outer', on='gene_id')

out_name = 'annotated_raw_counts.tsv'

merged_data.to_csv(out_name, sep='\t', index=False)
