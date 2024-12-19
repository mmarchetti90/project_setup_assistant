#!/usr/bin/env python3

"""
This script parses a filtered ATAC file and adds biotype
"""

### ---------------------------------------- ###

def parse_args():
    
    atac_file = argv[argv.index('--atac') + 1]

    gtf_file = argv[argv.index('--gtf') + 1]
    
    if '--chromosomes' in argv:
        
        chromosomes_of_interest = argv[argv.index('--chromosomes') + 1].split(',')
    
    else:
        
        chromosomes_of_interest = []
    
    return atac_file, gtf_file, chromosomes_of_interest

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

def add_info(atac_data, gtf_data):
    
    original_ids, new_ids, symbols, chromosomes, biotypes = [], [], [], [], []
    
    for gene in atac_data.index.values:
        
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
    
    #atac_columns = atac_data.columns.to_list()
    atac_columns = ['PeakID', 'Start', 'End', 'Length', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    
    atac_data = atac_data.assign(original_id = original_ids,
                                 gene_id = new_ids,
                                 gene_symbol = symbols,
                                 chromosome = chromosomes,
                                 biotype = biotypes)
    
    atac_data = atac_data.loc[:, ['original_id', 'gene_id', 'gene_symbol', 'chromosome', 'biotype'] + atac_columns]
    atac_data.index = range(atac_data.shape[0])
    atac_data.columns = ['original_id', 'gene_id', 'gene_symbol', 'chromosome', 'biotype', 'peak_id', 'peak_start', 'peak_end', 'peak_length', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    
    return atac_data

### ------------------MAIN------------------ ###

import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from sys import argv

### Parse CLI args

atac_file, gtf_file, chromosomes_of_interest = parse_args()

### Load ATAC data

atac = pd.read_csv(atac_file, sep='\t', index_col=0)

### Load GTF file

gtf = parse_gtf(gtf_file, chromosomes_of_interest)

### Add info to DEA file

atac = add_info(atac.copy(), gtf)

### Write to file

out_name = f'{atac_file.split("/")[-1][:-4]}_annotated.tsv'
atac.to_csv(out_name, sep='\t', index=False, header=True)

### Summarize DEA protein coding genes by chromosome

p_thr = 0.01
bad_biotypes = ['Mt_rRNA', 'Mt_tRNA', 'TEC', 'artifact', 'lncRNA',
                'miRNA', 'misc_RNA','rRNA', 'rRNA_pseudogene', 'ribozyme',
                'rRNA', 'rRNA_pseudogene', 'ribozyme', 'sRNA', 'scRNA',
                'scaRNA', 'snRNA', 'snoRNA', 'vault_RNA']
atac = atac.loc[(atac.padj < p_thr) &
                (~ atac.biotype.isin(bad_biotypes)),]

genes_by_chromosome = atac.groupby(by='chromosome').size()
if len(chromosomes_of_interest):
    
    genes_by_chromosome = genes_by_chromosome.loc[[chrom for chrom in chromosomes_of_interest if chrom in genes_by_chromosome.index]]

genes_by_chromosome = pd.DataFrame({'chromosome' : genes_by_chromosome.index,
                                    'gene_count' : genes_by_chromosome.values})

out_name = f'{atac_file.split("/")[-1][:-4]}_annotated_summary.tsv'
genes_by_chromosome.to_csv(out_name, sep='\t', index=False, header=True)

### Plot summary

out_name = f'{atac_file.split("/")[-1][:-4]}_protein_coding.png'
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