#!/usr/bin/env python3"""Averages the ATAC coverage profile around the TSSs of DEGs (all, upreg only, or downreg only) and plots themThe atac_manifest has the following shape (tab-separated):(N.B. samples and references have to be the same used for differential expression analysis)comparison_name    <output_prefix>samples            atac_sample_coverage_1.bed.gz,...,atac_sample_coverage_N.bed.gzreferences         atac_reference_coverage_1.bed.gz,...,atac_reference_coverage_N.bed.gz"""### ---------------------------------------- ###def parse_args():        # Parse manifest    atac_manifest_file = argv[argv.index('--atac_manifest') + 1]        # DEA file    dea_file = argv[argv.index('--dea') + 1]        # GTF file    gtf_file = argv[argv.index('--gtf') + 1]    # p value threshold for filtering DEGs    if '--p_thr' in argv:                p_thr = float(argv[argv.index('--p_thr') + 1])        else:                p_thr = 0.05        # log2FC threshold for filtering DEGs    if '--fc_thr' in argv:                fc_thr = float(argv[argv.index('--fc_thr') + 1])        else:                fc_thr = 0.        return atac_manifest_file, dea_file, gtf_file, p_thr, fc_thr### ---------------------------------------- ###def get_genes_coverage(file, targets, coords, ext=2500):        # Init coverage    coverage = np.zeros(ext * 2)        # Load sample coverage    col_dtypes = {0 : str, 1 : int, 2 : int, 3 : float}    cov = pd.read_csv(file, sep='\t', header=None, dtype=col_dtypes)    cov.columns = ['chrom', 'start', 'end', 'coverage']        # Subset coords for genes of interest    coords_sub = coords.loc[coords.gene_id.isin(targets),]        # Process genes    for _,(chrom, tss, strand, gene) in coords_sub.iterrows():            # Filter for coords        cov = cov.loc[(cov.chrom == str(chrom)) &                      (cov.start >= tss - ext) &                      (cov.end <= tss + ext),]            # Get gene coverage        gene_cov = np.zeros(ext * 2)        for row_n, row in cov.iterrows():                        chrom, start, end, c = row            gene_cov[int(start) - tss : int(end) - tss + 1] += c                # Reverse axis if strand is '-'        if strand == '-':                        gene_cov = np.flip(gene_cov)                # Average per number of genes        coverage += gene_cov / coords_sub.shape[0]        return coverage        def get_gene_coverage(coords, files, ext=2500):        # Init coverage file    coverages = np.zeros((len(files), ext * 2))        for f_num,f in enumerate(files):                # Load sample coverage        col_dtypes = {0 : str, 1 : int, 2 : int, 3 : float}        cov = pd.read_csv(f, sep='\t', header=None, dtype=col_dtypes)        cov.columns = ['chrom', 'start', 'end', 'coverage']                # Filter for coords        cov = cov.loc[(cov.chrom == str(coords[0])) &                      (cov.start >= coords[1] - ext) &                      (cov.end <= coords[1] + ext),]                # Fill group_coverage        for row_n, row in cov.iterrows():                        chrom, start, end, c = row            coverages[f_num, int(start) - coords[1] : int(end) - coords[1] + 1] = c        # Reverse axis if strand is '-'    if coords[2] == '-':                coverages = np.flip(coverages, axis=1)        return coverages### ---------------------------------------- ###def plot_averaged_coverage_data(sample_coverage, reference_coverage, gene_group='All'):        figsize = (10, 4)    ext = len(sample_coverage) // 2    #ymin, ymax = 0, max(sample_coverage.max(), reference_coverage.max())    fig, axs = plt.subplots(nrows=2,                            ncols=1,                            figsize=figsize,                            sharex=True,                            sharey=True)    # Plot coverage for each group    group_colors = ['red', 'blue']    for n, (data, color) in enumerate(zip([sample_coverage, reference_coverage], group_colors)):                axs[n].fill_between(x=np.arange(0, len(data), 1),                            y1=np.zeros(len(data)),                            y2=data,                            color=color,                            lw=0)                axs[n].plot(np.arange(0, len(data), 1),                    data,                    c=color,                    lw=1,                    marker=None,                    alpha=0.5)                #axs[n].set_ylim(ymin, ymax)        axs[n].set_ylabel('ATAC signal (CPM)', loc='center')                if n == 1:                        plt.xticks([0, ext, len(sample_coverage)], [-ext, 'TSS', ext])        plt.tight_layout()    plt.savefig(f'{gene_group}_averaged_profile.png', dpi=600)    plt.close()    ### ---------------------------------------- ###def plot_averaged_coverage_data_overlaid(sample_coverage, reference_coverage, gene_group='All'):        figsize = (10, 2)    ext = len(sample_coverage) // 2    #ymin, ymax = 0, max(sample_coverage.max(), reference_coverage.max())    plt.figure(figsize=figsize)        # Plot coverage for each group    group_colors = ['red', 'blue']    for n, (data, color) in enumerate(zip([sample_coverage, reference_coverage], group_colors)):                plt.fill_between(x=np.arange(0, len(data), 1),                         y1=np.zeros(len(data)),                         y2=data,                         color=color,                         alpha=0.5,                         lw=0)            #plt.set_ylim(ymin, ymax)    plt.ylabel('ATAC signal (CPM)', loc='center')    plt.xticks([0, ext, len(sample_coverage)], [-ext, 'TSS', ext])        plt.tight_layout()    plt.savefig(f'{gene_group}_averaged_profile_overlaid.png', dpi=600)    plt.close()### ------------------MAIN------------------ ###import numpy as npimport pandas as pdfrom matplotlib import pyplot as pltfrom sys import argv### Load data# Parse argsatac_manifest_file, dea_file, gtf_file, p_thr, fc_thr = parse_args()# Load ATAC manifestatac_files = {row.split('\t')[0] : row.split('\t')[1] if 'comparison_name' in row else row.split('\t')[1].split(',')              for row in open(atac_manifest_file, 'r').read().split('\n')              if len(row)}# DEA datadea = pd.read_csv(dea_file, sep='\t')dea.drop('original_id', axis=1, inplace=True)# Load GTFgtf = pd.read_csv(gtf_file, sep='\t', header=None, comment='#')gtf.columns = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']gtf = gtf.loc[gtf.feature == 'start_codon',]gene_ids = [[a.replace('gene_id ', '').replace('"', '').replace(';', '') for a in attr.split(';') if 'gene_id "' in a][0]            for attr in gtf.attributes.values]gtf = gtf.assign(gene_id = gene_ids)### Get list of DEGs and subset GTFdea = dea.loc[(dea.biotype == 'protein_coding') &              (dea.padj < p_thr) &              (dea.log2FoldChange.abs() > fc_thr),]gtf = gtf.loc[gtf.gene_id.isin(dea.gene_id.values), ['chrom', 'start', 'strand', 'gene_id']]### Create set of DEGsgene_sets = {'All' : dea.gene_id.values,              'Upregulated' : dea.loc[dea.log2FoldChange > 0, 'gene_id'].values,             'Downregulated' : dea.loc[dea.log2FoldChange < 0, 'gene_id'].values}### Extract coveragesfor set_name, genes in gene_sets.items():        # Get coverage for samples    sample_coverage = []    for sample in atac_files['samples']:                sample_coverage.append(get_genes_coverage(sample, genes, gtf))        # Average across samples    sample_coverage = np.mean(sample_coverage, axis=0)    # Smooth    sample_coverage = np.array([sample_coverage[max(0, i - 10) : i + 10].mean() for i in range(sample_coverage.shape[0])])        # Get coverage for controls    reference_coverage = []    for ref in atac_files['references']:                reference_coverage.append(get_genes_coverage(ref, genes, gtf))        # Average across samples    reference_coverage = np.mean(reference_coverage, axis=0)        # Smooth    reference_coverage = np.array([reference_coverage[max(0, i - 10) : i + 10].mean() for i in range(reference_coverage.shape[0])])        # Plot data    plot_averaged_coverage_data(sample_coverage, reference_coverage, gene_group=set_name)    plot_averaged_coverage_data_overlaid(sample_coverage, reference_coverage, gene_group=set_name)    