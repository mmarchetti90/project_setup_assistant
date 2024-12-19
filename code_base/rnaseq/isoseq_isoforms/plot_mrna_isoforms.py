#!/usr/bin/env python3

"""
This script plots isoforms for individual genes
Supports multi-sample pigeon classification files
"""

### ---------------------------------------- ###

def pars_argv():
    
    # Gene names provided as a comma separated list, e.g. Egfr,Egf,Wnt6,Col4a1
    genes = argv[argv.index('--gene') + 1].split(',')
    
    # Path to Isoseq GFF file from prepared by Pigeon
    gff_path = argv[argv.index('--gff') + 1]
    gff = pd.read_csv(gff_path, sep='\t', comment='#', header=None)
    gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    
    # Path to Pigeon classification file
    classification_path = argv[argv.index('--classification') + 1]
    classification = pd.read_csv(classification_path, sep='\t')
    
    return genes, gff, classification

### ---------------------------------------- ###

def plot_gene(gene, gff, classification):
    
    if gene not in classification.associated_gene.values:
        
        print(f'ERROR: gene {gene} not found.')
        return ''
    
    # Extracting isoforms
    isoforms_ids = classification.loc[classification.associated_gene == gene, 'isoform'].values
    gene_pb_id = '.'.join(isoforms_ids[0].split('.')[:-1])
    
    # Extracting exons
    attributes = {f'gene_id "{".".join(iso.split(".")[:-1])}"; transcript_id "{iso}";' : iso for iso in isoforms_ids}
    #attributes = {f'gene_id "{gene_pb_id}"; transcript_id "{iso}";' : iso for iso in isoforms_ids}
    gene_exons = gff.loc[(gff.attribute.isin(attributes.keys())) & (gff.feature == 'exon'),].copy()
    gene_exons['isoform'] = [attributes[iso] for iso in gene_exons.attribute]
    
    # Get count columns and sample names
    counts_cols = [col for col in classification.columns if col.startswith('FL_TPM') and not col.endswith('_log10')]
    samples = [col.replace('FL_TPM.', '') for col in counts_cols]
    if len(samples) == 1 and samples[0] == '':
        
        samples = ['Sample']
    
    # Init plot and subplots
    figsize = (15 + 1 * len(samples), max(3, 1 * len(isoforms_ids) + 2))
    f, (isoform_plots, expression_plot) = plt.subplots(1, 2, figsize=figsize, width_ratios=[15, 1 * len(samples)])
    
    # Plot isoforms
    plt.subplot(1, 2, 1)
    exon_linewidth, exon_color = 14, 'red'
    intron_linewidth, intron_color = 2, 'black'
    ymin, ymax = 0, len(isoforms_ids) + 1
    xlab = gene_exons.seqname.values[0]
    ylab = f'{gene}\nisoforms'
    for iso_n,iso in enumerate(isoforms_ids):
        
        # Define y coordinate
        y = len(isoforms_ids) - iso_n
        
        # Extract isoform exons and introns
        iso_exons = gene_exons.loc[gene_exons.isoform == iso, ['start', 'end']].values.tolist()
        iso_exons.sort(key=lambda i: i[0])
        iso_introns = [[ex1[1], ex2[0]] for ex1,ex2 in zip(iso_exons[:-1], iso_exons[1:])]
        
        # Plot introns
        for intron in iso_introns:
            
            isoform_plots.hlines(y, intron[0], intron[1], linewidth=intron_linewidth, colors=intron_color)
    
        # Plot exons
        for exon in iso_exons:
            
            isoform_plots.hlines(y, exon[0], exon[1], linewidth=exon_linewidth, colors=exon_color)
    
    plt.xlabel(xlab, fontweight='bold')
    plt.ylabel(ylab, fontweight='bold')
    plt.ylim((ymin, ymax))
    plt.yticks(ticks=range(len(isoforms_ids), 0, - 1), labels=isoforms_ids, fontweight='bold')
    #isoform_plots.spines.top.set_visible(False)
    #isoform_plots.spines.right.set_visible(False)
    #isoform_plots.spines.left.set_visible(False)
    
    # Plot isoforms' expression
    """
    # Seaborn heatmap
    plt.subplot(1, 2, 2)
    heatmap_data = classification.loc[classification.associated_gene == gene, counts_cols].copy()
    heatmap_data.index = isoforms_ids
    heatmap_data.columns = samples
    sns.heatmap(data=heatmap_data,
                cmap='viridis',
                linewidth=0.5,
                linecolor='black',
                square=True,
                xticklabels=True,
                yticklabels=False)
    plt.xticks(fontweight='bold', rotation=45, ha='right')
    """
    # Dot plot
    plt.subplot(1, 2, 2)
    xmin, xmax = 0, len(samples) + 1
    ymin, ymax = 0, len(isoforms_ids) + 1
    x_coords = [i + 1 for _ in range(len(isoforms_ids)) for i in range(len(samples))]
    y_coords = [i for i in range(len(isoforms_ids), 0, - 1) for _ in range(len(samples))]
    expression = classification.loc[classification.associated_gene == gene, counts_cols].values.flatten()
    expression_plot_fig = expression_plot.scatter(x=x_coords, y=y_coords, s=500, c=expression, marker='o', cmap='viridis', edgecolors='black', linewidths=1)
    plt.title('Expression (TPM)', fontweight='bold')
    plt.xlabel('')
    plt.ylabel('')
    plt.xlim((xmin, xmax))
    plt.ylim((ymin, ymax))
    plt.xticks(ticks=range(1, len(samples) + 1, 1), labels=samples, fontweight='bold', rotation=90)
    plt.yticks([], [])
    plt.box(False)
    
    f.colorbar(expression_plot_fig, ax=expression_plot, orientation='vertical', fraction=0.05)
    
    #plt.subplots_adjust(wspace=0.1)
    plt.tight_layout()
    plt.savefig(f'{gene}_isoforms.png', dpi=300)
    plt.close()
    
### ------------------MAIN------------------ ###

import pandas as pd
#import seaborn as sns

from matplotlib import pyplot as plt
from sys import argv

### Parse cli

genes, gff, classification = pars_argv()

### Plot individual genes

for gene in genes:
    
    plot_gene(gene, gff, classification)
