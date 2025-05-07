#!/usr/bin/env python3

"""
This script parses a DEA file from run_deseq2.R or run_dexseq.R and adds useful info
"""

### ---------------------------------------- ###

def parse_args():
    
    dea_file = argv[argv.index('--dea') + 1]

    gtf_file = argv[argv.index('--gtf') + 1]
    
    if '--level' not in argv:
        
        analysis_level = 'gene'
    
    else:
        
        analysis_level = argv[argv.index('--level') + 1]
    
    if '--chromosomes' in argv:
        
        chromosomes_of_interest = argv[argv.index('--chromosomes') + 1].split(',')
    
    else:
        
        chromosomes_of_interest = []
    
    if '--biotypes' in argv:
        
        biotypes_of_interest = argv[argv.index('--biotypes') + 1].split(',')
    
    else:
        
        biotypes_of_interest = []
    
    return dea_file, gtf_file, analysis_level, chromosomes_of_interest, biotypes_of_interest

### ---------------------------------------- ###

def parse_gtf(path, desired_biotypes=[], desired_chromosomes=[], analysis_level='gene'):
    
    # Load GTF
    
    gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    gtf_data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Only keep genes and exons

    gtf_data = gtf_data.loc[gtf_data.feature.isin(['gene', 'exon']), ['seqname', 'start', 'end', 'strand', 'attribute']]
    
    # Get biotype and gene id

    biotypes, gene_ids, gene_symbols, transcript_ids, exon_ids = [], [], [], [], []
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
        
        if 'transcript_id' in info:
            
            transcript = info[info.index('transcript_id "') + len('transcript_id "'):]
            transcript = transcript[:transcript.index('"')]
        
        else:
            
            transcript = ''
        
        transcript_ids.append(transcript)
        
        if 'exon_id' in info:
            
            exon = info[info.index('exon_id "') + len('exon_id "'):]
            exon = exon[:exon.index('"')]
        
        else:
            
            exon = ''
        
        exon_ids.append(exon)

    gtf_data['biotype'] = biotypes
    gtf_data['gene_id'] = gene_ids
    gtf_data['gene_symbol'] = gene_symbols
    gtf_data['transcript_id'] = transcript_ids
    gtf_data['exon_id'] = exon_ids
    
    # Fix dtypes
    
    gtf_data[['start', 'end']] = gtf_data[['start', 'end']].astype(int)
    
    # Filter based on biotype
    
    if len(desired_biotypes):

        gtf_data = gtf_data.loc[gtf_data.biotype.isin(desired_biotypes),]

    # Filter for desired chromosomes

    if len(desired_chromosomes):

        gtf_data = gtf_data.loc[gtf_data.seqname.isin(desired_chromosomes),]
    
    # Remove genes without gene_symbol
    
    #gtf_data = gtf_data.loc[gtf_data['gene_symbol'] != '',]
    
    # Remove undesired levels
    
    if analysis_level == 'gene':
        
        gtf_data.drop(['attribute', 'transcript_id', 'exon_id'], axis='columns', inplace=True)
        gtf_data.drop_duplicates(inplace=True)
    
    elif analysis_level == 'transcript':
        
        gtf_data.drop(['attribute', 'exon_id'], axis='columns', inplace=True)
        gtf_data.drop_duplicates(inplace=True)
    
    else:
        
        gtf_data.drop(['attribute'], axis='columns', inplace=True)
        gtf_data.drop_duplicates(inplace=True)
    
    return gtf_data

### ---------------------------------------- ###

def add_info(dea_data, gtf_data, analysis_level='gene'):
    
    annotation_header = ['original_id', 'gene_id', 'gene_symbol', 'transcript_id', 'exon_id', 'chromosome', 'element_start', 'element_end', 'strand', 'biotype']
    annotation = []
    
    for original_id in dea.index.values:
        
        if analysis_level == 'gene':
            
            gtf_data_sub = gtf_data.loc[(gtf_data.gene_id == original_id.split('.')[0]) | (gtf_data.gene_symbol == original_id.split('.')[0]),]
            
        elif analysis_level == 'transcript':
            
            gtf_data_sub = gtf_data.loc[gtf_data.transcript_id == original_id,]
        
        else:
            
            gtf_data_sub = gtf_data.loc[gtf_data.exon_id == original_id,]
        
        # Skip if element is not found in the filtered gtf (i.e. element is not in desired chromosomes or has an unwanted biotype)
        if gtf_data_sub.shape[0] == 0:
            
            continue
        
        gene_id = gtf_data_sub['gene_id'].values[0]
        gene_symbol = gtf_data_sub['gene_symbol'].values[0]
        transcript_id = '' if analysis_level in ['gene'] else gtf_data_sub['transcript_id'].values[0]
        exon_id = '' if analysis_level in ['gene', 'transcript'] else gtf_data_sub['transcript_id'].values[0]
        chrom = gtf_data_sub['seqname'].values[0]
        start = gtf_data_sub['start'].values.min()
        end = gtf_data_sub['end'].values.max()
        strand = gtf_data_sub['end'].values[0]
        biotype = gtf_data_sub['biotype'].values[0]
        
        annotation.append([original_id, gene_id, gene_symbol, transcript_id, exon_id, chrom, start, end, strand, biotype])
        
    annotation = pd.DataFrame(annotation, columns=annotation_header)
    
    if analysis_level == 'gene':
        
        annotation.drop(['transcript_id', 'exon_id'], axis='columns', inplace=True)
        annotation.drop_duplicates(inplace=True)
    
    elif analysis_level == 'transcript':
        
        annotation.drop(['exon_id'], axis='columns', inplace=True)
        annotation.drop_duplicates(inplace=True)
    
    else:
        
        annotation.drop_duplicates(inplace=True)
    
    dea_data = pd.merge(annotation, dea_data, how='inner', left_on='original_id', right_index=True)
    dea_data.drop(['original_id'], axis='columns', inplace=True)
    
    return dea_data

### ------------------MAIN------------------ ###

import pandas as pd
import re
import seaborn as sns

from matplotlib import pyplot as plt
from sys import argv

### Parse CLI args

dea_file, gtf_file, analysis_level, chromosomes_of_interest, biotypes_of_interest = parse_args()

### Load data

# DEA data

if analysis_level != 'exon':
    
    dea = pd.read_csv(dea_file, sep='\t', index_col=0)

else:
    
    dea = pd.read_csv(dea_file, sep='\t', index_col=1)

# Load GTF file

gtf = parse_gtf(gtf_file, biotypes_of_interest, chromosomes_of_interest, analysis_level)

### Add info to DEA file

dea = add_info(dea.copy(), gtf, analysis_level)

### Write to file

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated.tsv'
dea.to_csv(out_name, sep='\t', index=False, header=True)

### Summarize genes by chromosome

chromosomes_of_interest = list(set(dea.chromosome))
chromosomes_of_interest.sort()

elements_by_chromosome = dea.groupby(by='chromosome').size()
elements_by_chromosome = elements_by_chromosome.loc[[chrom for chrom in chromosomes_of_interest if chrom in elements_by_chromosome.index]]
elements_by_chromosome = pd.DataFrame({'chromosome' : elements_by_chromosome.index,
                                       f'{analysis_level}_count' : elements_by_chromosome.values})

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated_summary.tsv'
elements_by_chromosome.to_csv(out_name, sep='\t', index=False, header=True)

### Plot summary

log2fc_col = [col for col in dea.columns.values if col.lower().startswith('log2f')][0] # Log2FC columns are named log2FoldChange from DESeq2, and log2fold_<comparison> for DEXSeq

dea.sort_values(by='chromosome', inplace=True)

# All genes

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated_summary_allgenes.png'
figsize = (len(set(elements_by_chromosome.chromosome)), 5)
plt.figure(figsize=figsize)
strip = sns.stripplot(dea, x='chromosome', y=log2fc_col, hue='biotype', edgecolor='black', linewidth=0.25)
plt.xlabel('Chromosome')
plt.ylabel('$log_2FC$')
plt.xticks(rotation=45, ha='right')
sns.move_legend(strip, loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig(out_name, dpi=300)
plt.close()

# p < 0.05

out_name = f'{dea_file.split("/")[-1][:-4]}_annotated_summary_p-0.05.png'
figsize = (len(set(elements_by_chromosome.chromosome)), 5)
plt.figure(figsize=figsize)
strip = sns.stripplot(dea.loc[dea.padj < 0.05,], x='chromosome', y=log2fc_col, hue='biotype', edgecolor='black', linewidth=0.25)
plt.xlabel('Chromosome')
plt.ylabel('$log_2FC$')
plt.xticks(rotation=45, ha='right')
sns.move_legend(strip, loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig(out_name, dpi=300)
plt.close()
