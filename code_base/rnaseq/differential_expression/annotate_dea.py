#!/usr/bin/env python3

"""
This script parses a DEA file from run_deseq2.R and adds useful info
"""

### ---------------------------------------- ###

def parse_args():
    
    dea_file = argv[argv.index('--dea') + 1]

    gtf_file = argv[argv.index('--gtf') + 1]
    
    if '--chromosomes' in argv:
        
        chromosomes_of_interest = argv[argv.index('--chromosomes') + 1].split(',')
    
    else:
        
        chromosomes_of_interest = []
    
    return dea_file, gtf_file, chromosomes_of_interest

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

### ---------------------------------------- ###

def add_info(dea_data, gtf_data):
    
    original_ids, new_ids, symbols, chromosomes, biotypes = [], [], [], [], []
    
    for gene in dea.index.values:
        
        original_ids.append(gene)
        
        gene = gene.split('.')[0]
        
        if gene in gtf_data.gene_id.values:
            
            gene_id = gene
            
            _, gene_symbol, gene_chromosome, gene_biotype = gtf_data.loc[gtf_data.gene_id == gene,].iloc[0,].values
        
        elif gene in gtf_data.gene_symbol.values:
            
            gene_symbol = gene
            
            gene_id, _, gene_chromosome, gene_biotype = gtf_data.loc[gtf_data.gene_symbol == gene,].iloc[0,].values
        
        else:
            
            gene_id, gene_symbol, gene_chromosome, gene_biotype = ['Unknown' for _ in range(4)]
        
        new_ids.append(gene_id)
        symbols.append(gene_symbol)
        chromosomes.append(gene_chromosome)
        biotypes.append(gene_biotype)
    
    dea_columns = dea_data.columns.to_list()
    
    dea_data = dea_data.assign(original_id = original_ids,
                               gene_id = new_ids,
                               gene_symbol = symbols,
                               chromosome = chromosomes,
                               biotype = biotypes)
    
    dea_data = dea_data.loc[:, ['original_id', 'gene_id', 'gene_symbol', 'chromosome', 'biotype'] + dea_columns]
    dea_data.index = range(dea_data.shape[0])
    
    return dea_data

### ------------------MAIN------------------ ###

import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from sys import argv

### Parse CLI args

dea_file, gtf_file, chromosomes_of_interest = parse_args()

### Load DEA data

dea = pd.read_csv(dea_file, sep='\t', index_col=0)

### Load GTF file

gtf = parse_gtf(gtf_file, chromosomes_of_interest)

### Add info to DEA file

dea = add_info(dea.copy(), gtf)

### Write to file

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated.tsv'
dea.to_csv(out_name, sep='\t', index=False, header=True)

### Summarize genes by chromosome

genes_by_chromosome = dea.groupby(by='chromosome').size()
genes_by_chromosome = genes_by_chromosome.loc[[chrom for chrom in chromosomes_of_interest if chrom in genes_by_chromosome.index]]
genes_by_chromosome = pd.DataFrame({'chromosome' : genes_by_chromosome.index,
                                    'gene_count' : genes_by_chromosome.values})

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated_summary.tsv'
genes_by_chromosome.to_csv(out_name, sep='\t', index=False, header=True)

### Plot summary

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated_summary.png'
figsize = (len(set(genes_by_chromosome.chromosome)) / 2, 2.5)
plt.figure(figsize=figsize)
sns.barplot(genes_by_chromosome,
            x='chromosome',
            y='gene_count',
            lw=1,
            edgecolor='black')
plt.xlabel('Chromosome')
plt.ylabel('Gene count')
plt.tight_layout()
plt.savefig(out_name, dpi=300)
plt.close()