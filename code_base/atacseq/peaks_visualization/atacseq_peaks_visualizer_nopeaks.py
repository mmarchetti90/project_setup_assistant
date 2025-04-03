#!/usr/bin/env python3

"""
This script will plot peaks from atacseq
Input is a data manifest, a tsv file with following lines:
- "annotation	</path/to/annotation/gtf" to declare GTF file for annotation;
- "sample <sample_group> </path/to/replicate/of/sample/group/bed>" to declare a bed file belonging to a specific sample group
- "gene <gene_id or gene_name>" to declare a gene to be plotted
N.B. Samples belonging to the same group will be averaged.
"""

### ---------------------------------------- ###

def parse_args():
    
    # Parse manifest
    data_manifest_file = argv[argv.index('--data_manifest') + 1]
    
    samples_path = {}
    genes_list = []
    
    with open(data_manifest_file) as info:
        
        for line in info:
            
            info = line.replace('\n', '').split('\t')
            
            if info[0] == 'annotation':
                
                annotation_path = info[1]
            
            elif info[0] == 'sample':
                
                try:
                    
                    samples_path[info[1]].append(info[2])
                    
                except:
                    
                    samples_path[info[1]] = [info[2]]
            
            elif info[0] == 'gene':
                
                genes_list.append(info[1])
            
            else:
                
                pass
    
    # Genes region extension
    if '--gene_ext' in argv:
        
        gene_extension = int(argv[argv.index('--gene_ext') + 1])
    
    else:
        
        gene_extension = 1000
    
    return annotation_path, samples_path, genes_list, gene_extension

### ---------------------------------------- ###

def get_gene_structure(ann, gn):
    
    # Subset for gene of interest
    ann = ann.loc[[f'gene_id "{gn}"' in attr or f'gene_name "{gn}"' in attr for attr in ann.attributes.values], ]
    
    # Get gene coords
    coords = ann.loc[ann.feature == 'gene', ['chrom', 'start', 'end']].values[0]
    
    # Keep only info for gene structure (exons and UTRs)
    ann = ann.loc[ann.feature.isin(['five_prime_utr', 'three_prime_utr', 'exon'])]
    
    return coords, ann

### ---------------------------------------- ###

def get_gene_coverage(coords, s_p):
    
    coverages = {}
    for group, samples in s_p.items():
        
        group_coverage = np.zeros((len(samples), coords[2] - coords[1] + 1))

        for s_num,sample in enumerate(samples):
            
            # Load sample coverage
            col_dtypes = {0 : str, 1 : int, 2 : int, 3 : float}
            sample_coverage = pd.read_csv(sample, sep='\t', header=None, dtype=col_dtypes)
            sample_coverage.columns = ['chrom', 'start', 'end', 'coverage']
            
            # Filter for coords
            sample_coverage = sample_coverage.loc[(sample_coverage.chrom == str(coords[0])) &
                                                  (sample_coverage.start >= coords[1]) &
                                                  (sample_coverage.end <= coords[2]),]
            
            sample_coverage = sample_coverage.loc[(sample_coverage.chrom == str(coords[0])) &
                                                  (sample_coverage.start <= coords[2]) &
                                                  (sample_coverage.end >= coords[1]),]
            
            # Fill group_coverage
            for row_n, row in sample_coverage.iterrows():
                
                chrom, start, end, cov = row
                group_coverage[s_num, int(start) - coords[1] : int(end) - coords[1] + 1] = cov

        coverages[group] = group_coverage.copy()
    
    return coverages

### ---------------------------------------- ###

def plot_individual_coverage_data(gn, coords, struct, covs):
    
    figsize = (10, sum([len(c) for c in covs.values()]) + 1)
    height_ratios = [2 for _ in range(sum([len(c) for c in covs.values()]))] + [0.55]
    xmin, xmax = coords[1:]
    ymin, ymax = 0, max([c.max() for c in covs.values()])

    fig, axs = plt.subplots(nrows=sum([len(c) for c in covs.values()]) + 1,
                            ncols=1,
                            figsize=figsize,
                            height_ratios=height_ratios,
                            sharex=True)

    plt.xlim(xmin, xmax)

    # Plot coverage for each group
    group_colors = [colors.to_hex(color) for color in plt.cm.tab20.colors]
    
    n = 0

    for n_group, ((group, group_cov), color) in enumerate(zip(covs.items(), group_colors)):
        
        for n_sample, sample_cov in enumerate(group_cov):
        
            axs[n].fill_between(x=np.arange(coords[1], coords[2] + 1, 1),
                                y1=np.zeros(coords[2] - coords[1] + 1),
                                y2=sample_cov,
                                color=color,
                                lw=0)
            
            axs[n].plot(np.arange(coords[1], coords[2] + 1, 1),
                        sample_cov,
                        c=color,
                        lw=1,
                        marker=None,
                        alpha=0.5)
            
            axs[n].set_ylim(ymin, ymax)
            axs[n].set_title(f'{group} {n_sample+1}', loc='center')
            axs[n].axis('off')
            
            n += 1

    # Plot gene body
    intron_size = 1

    axs[n].plot([struct.start.min(), struct.end.max() + 1],
                [0, 0],
                c='black',
                lw=intron_size)

    # Plot gene elements
    utr_size = 5
    exon_size = 8

    for _,element in struct.iterrows():
        
        if element.feature == 'exon':
            
            axs[n].plot([element.start, element.end + 1],
                        [0, 0],
                        c='black',
                        lw=exon_size)

        else:
            
            axs[n].plot([element.start, element.end + 1],
                        [0, 0],
                        c='black',
                        lw=utr_size)

    # Plot arrows indicating direction of gene
    strand = struct.strand.values[0]

    arrows_x = np.linspace(struct.start.min(), struct.end.max() + 1, 30)
    if strand == '-' :
        
        arrows_x = arrows_x[::-1]
        
    #axs[n].plot(arrows_x,
    #            np.zeros(arrows_x.shape[0]),
    #            c='black',
    #            lw=0,
    #            marker='>' if strand == '+' else '<',
    #            markersize=3)

    #axs[n].plot(arrows_x,
    #            np.zeros(arrows_x.shape[0]),
    #            c='white',
    #            lw=0,
    #            marker='>' if strand == '+' else '<',
    #            markersize=2)

    axs[n].plot(arrows_x,
                np.zeros(arrows_x.shape[0]),
                c='black',
                lw=0,
                marker='4' if strand == '+' else '3',
                markersize=8)

    axs[n].axis('off')

    plt.tight_layout()
    plt.savefig(f'{gene}_profile.png', dpi=600)
    plt.close()
    
### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from matplotlib import colors
from matplotlib import pyplot as plt
from sys import argv

### Parse data manifest

annotation_path, samples_path, genes_list, gene_extension = parse_args()

### Load gtf file

gtf = pd.read_csv(annotation_path, sep='\t', header=None, comment='#')
gtf.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

### Plot genes

for gene in genes_list:
    
    # Extract gene structure
    gene_coords, gene_structure = get_gene_structure(gtf, gene)
    
    # Extend region
    gene_coords[1] -= gene_extension
    gene_coords[2] += gene_extension
    
    # Extract coverage in gene region
    gene_coverages = get_gene_coverage(gene_coords, samples_path)
    
    # Plot
    plot_individual_coverage_data(gene, gene_coords, gene_structure, gene_coverages.copy())
