#!/usr/bin/env python3

"""
This script parses a differential analysis (e.g. expression, accessibility, etc) and looks for
over-representation or depletion of terms (e.g. GO terms)
"""

### ---------------------------------------- ###

def print_help():
    
    """
    This function displays a description of the script
    """
    
    help_text = """
    This script parses a differential analysis (e.g. expression, accessibility, etc) and looks for
    over-representation or depletion of terms (e.g. GO terms)

    Parameters
    ----------
    --help : optional
        Prints this help text and quits.
    --differential_analysis : string path, required
        Path to differential analysis.
        Must be a tab-delimited file with 3 columns: gene_symbol, log2FoldChange, and padj
    --term2gene : string path, required
        Path to the terms annotations.
        Must be a tab-delimited file with 2 columns: category, gene_symbol
    --universe : string path, optional
        Path to a txt list of genes to be used as background.
        If not specified, all the genes in term2gene will be used.
    --padj : float, optional
        Threshold for differential_analysis padj.
        Default = 0.05
    --log2fc : float, optional
        Threshold for differential_analysis log2FC.
        Default = 0
    --min_genes : int, optional
        Minumum number of genes in a category.
        Default = 10
    --max_genes : int, optional
        Maximum number of genes in a category.
        Default = 500
    --filter_genes : bool, optional
        If True, the differential_analysis is filtered to remove genes not in the universe (like in clusterProfiler).
        If False, genes in differential_analysis missing from the universe are added to the latter.
        Added for compatibility.
        Default = False
    --output_prefix : string, optional
        Prefix for output files.
        Default = 'enrichment''
    --split_upreg_downreg : bool, optional
        If True, enrichment will be tested using 3 datasets: all genes, only upregulated ones, and only downregulated ones.
        If False, enrichment will be tested with all genes only.
        Default = True
    --enrichment_only : optional
       If specified, only enrichemnt will be tested.
       Mutually exclusive with --depletion_only.
    --depletion_only : optional
       If specified, only depletion will be tested.
       Mutually exclusive with --enrichment_only.
    """

    print(help_text)

### ---------------------------------------- ###

def parse_args():
    
    print(argv)
    
    # p-value threshold
    if '--padj' in argv:
        
        p_thr = float(argv[argv.index('--padj') + 1])
    
    else:
        
        p_thr = 0.05
    
    # log2FC threshold
    if '--log2fc' in argv:
        
        log2fc_thr = float(argv[argv.index('--log2fc') + 1])
    
    else:
        
        log2fc_thr = 0.
        
    # Load differential analysis
    da_path = argv[argv.index('--differential_analysis') + 1]
    da = pd.read_csv(da_path, sep='\t')
    da = da.loc[(da.log2FoldChange.abs().values > log2fc_thr) & (da.padj < p_thr),]
    da = da.loc[~da.gene_symbol.isna(),]
    
    # Load term2gene
    t2g_path = argv[argv.index('--term2gene') + 1]
    t2g = pd.read_csv(t2g_path, sep='\t')
    t2g.drop_duplicates(inplace=True)
    t2g.columns = ['category', 'gene_symbol']
    
    if '--universe' in argv:
        
        u_path = argv[argv.index('--universe') + 1]
        u = open(u_path, 'r').read().split('\n')
    
    else:
        
        u = []

    # Minimum size of gene categories
    if '--min_genes' in argv:
        
        min_g = int(argv[argv.index('--min_genes') + 1])
    
    else:
        
        min_g = 10

    # Maximum size of gene categories
    if '--max_genes' in argv:
        
        max_g = int(argv[argv.index('--max_genes') + 1])
    
    else:
        
        max_g = 500

    # Remove genes in differential_analysis missing from the universe?
    if '--filter_genes' in argv:
        
        filt_g = (argv[argv.index('--filter_genes') + 1].lower() == "true")
    
    else:
        
        filt_g = False
    
    # Prefix of output files
    if '--output_prefix' in argv:
        
        out_p = argv[argv.index('--output_prefix') + 1]
    
    else:
        
        out_p = 'enrichment'
    
    if '--split_upreg_downreg' in argv:
        
        split_g = (argv[argv.index('--split_upreg_downreg') + 1].lower() == "true")
    
    else:
        
        split_g = True

    # Enrichment type
    if '--enrichment_only' in argv:

        e_type = 'enrichment'

    elif '--depletion_only' in argv:

        e_type = 'depletion'

    else:

        e_type = 'both'
    
    return da, t2g, u, min_g, max_g, filt_g, out_p, split_g, e_type

### ---------------------------------------- ###

class enrichment_analysis:

    """
    Class for the over-representation or depletion analysis of gene sets.

    Parameters
    ----------
    gene_list : list of strings
        Genes to be tested for enrichment.
    term2gene : pandas dataframe
        Dataframe with the following 2 columns:
        category = names of categories to test for enrichment/depletion
        gene_symbol = symbols of genes belonging to a specific category (1 gene per row)
    universe : list of strings, optional
        Genes to be used as background.
        If not provided, all genes in term2gene will be used.
        Note that genes in gene_list will be added if missing.
        Default = []
    output_prefix : string, optional
        Prefix for output files.
        Default = 'enrichment''
    filter_genes : bool, optional
        If True, genes in differential_analysis missing from the universe are removed.
        If False, they are added to the universe
        Default = False

    Attributes
    ----------
    results : pd.DataFrame
        Results of the analysis.
        Created by the process_data function

    Methods
    -------
    process_data()
        Processes the data and creates a results table.
    plot_enrichment()
        Plots the results in a bubble plot.
    """

    def __init__(self, gene_list, term2gene, universe=[], output_prefix='enrichment', filter_genes=False):

        # Set class vars
        self.output_prefix = output_prefix
        self.gene_list = np.array(gene_list)
        self.term2gene = term2gene

        # Define universe
        if filter_genes:

            if not len(universe):
            
                self.universe = np.unique(term2gene.gene_symbol.values)
                self.gene_list = self.gene_list[np.isin(self.gene_list, self.universe)]
                self.term2gene = self.term2gene.loc[self.term2gene.gene_symbol.isin(self.universe),]
            
            else:
                
                self.universe = np.unique(universe)
                self.gene_list = self.gene_list[np.isin(self.gene_list, self.universe)]

        else:
        
            if not len(universe):
                
                self.universe = np.unique(np.concatenate([term2gene.gene_symbol.values, gene_list]).astype(str))
            
            else:
                
                self.universe = np.unique(np.concatenate([universe, term2gene.gene_symbol.values, gene_list]).astype(str))

    ### ------------------------------------ ###
    ### PROCESSING                           ###
    ### ------------------------------------ ###
    
    def process_data(self, min_genes=10, max_genes=500, enrichment_type='both'):
        
        """
        This function runs the enrichment/depletion analysis.
        
        min_genes : int, optional
            Minumum number of genes in a category.
            Default = 5
        max_genes : int, optional
            Maximum number of genes in a category.
            Default = 500
        enrichment_type : string, optional
            Set to 'both' to compute both enrichment and depletion, or set to 'enrichment' or 'depletion', otherwise.
            Default = 'both'
        """
        
        # Filter term2gene
        category_sizes = pd.Series(self.term2gene.groupby(by='category').size())
        good_categories = category_sizes.index[(category_sizes >= min_genes) & (category_sizes <= max_genes)].values
        term2gene = self.term2gene.loc[self.term2gene.category.isin(good_categories)]
        
        # Init constants
        list_size, universe_size = len(self.gene_list), len(self.universe)
        
        # Init results table
        results = {'category' : [],
                   'gene_ratio' : [],
                   'background_ratio' : [],
                   'actual_hits' : [],
                   'expected_hits' : [],
                   'enrichment_score' : [],
                   'behaviour' : [],
                   'pval' : [],
                   'padj' : [],
                   'gene_hits' : []}
        
        # Test categories
        for n,category in enumerate(np.unique(term2gene.category.values)):
            
            if (n + 1) % 100 == 0 and n > 0:
                
                print(f'Processed {n + 1} / {len(good_categories)} categories', end='\r')
            
            # Get genes in the category in gene_list and the universe
            category_genes = np.unique(term2gene.gene_symbol[term2gene.category == category].values)
            list_hits = '/'.join(category_genes[np.isin(category_genes, self.gene_list)])
            list_hits_num = np.isin(category_genes, self.gene_list).sum()
            universe_hits_num = len(category_genes)
            
            # Write ratios and calculate expected number of genes
            gene_ratio, bg_ratio = f'{list_hits_num}/{list_size}', f'{universe_hits_num}/{universe_size}'
            expected_hits = round(list_size * (universe_hits_num / universe_size), 0)
            enrichment_score = np.log2((list_hits_num / list_size) / (universe_hits_num / universe_size))
            
            # Testing enrichment/depletion
            if (list_hits_num >= expected_hits and enrichment_type == 'both') or (enrichment_type == 'enrichment'):
                
                # Enrichment
                behaviour = 'enrichment'
                
                # Hypergeometric test
                # This was added for flavor, it's essentially the same as above
                #pval = self.hypergeom_testing(list_hits_num, list_size, universe_hits_num, universe_size)
                
                # Fisher's exact test
                pval = self.fisher_testing(list_hits_num, list_size - list_hits_num, universe_hits_num - list_hits_num, universe_size - list_size - universe_hits_num + list_hits_num)
            
            else:
                
                # Depletion
                behaviour = 'depletion'
                
                # Hypergeometric test
                # This was added for flavor, it's essentially the same as above
                #pval = self.hypergeom_testing(list_size - list_hits_num, list_size, universe_size - universe_hits_num, universe_size)
                
                # Fisher's exact test
                pval = self.fisher_testing(list_size - list_hits_num, list_hits_num, universe_size - list_size - universe_hits_num + list_hits_num, universe_hits_num - list_hits_num)
            
            results['category'].append(category)
            results['gene_ratio'].append(gene_ratio)
            results['background_ratio'].append(bg_ratio)
            results['actual_hits'].append(list_hits_num)
            results['expected_hits'].append(expected_hits)
            results['enrichment_score'].append(enrichment_score)
            results['behaviour'].append(behaviour)
            results['pval'].append(pval)
            results['gene_hits'].append(list_hits)
        
        # Adjustting p-value
        results['padj'] = fdrcorrection(results['pval'], alpha=0.05, is_sorted=False)[1]
        
        # Convert to dataframe and sort
        results = pd.DataFrame(results)
        results.sort_values(by='padj', ascending=True, inplace=True)
        
        # Save to file
        results.to_csv(f'{self.output_prefix}.tsv', sep='\t', index=False, header=True)
        
        self.results = results
    
    ### ------------------------------------ ###
    
    @staticmethod
    def hypergeom_testing(k, N, n, M):
        
        """
        Hypergeometric test for enrichment
        
        Parameters
        ----------
        k : int
            Number of significant genes in the category
        N : int
            Number of significant genes
        n : int
            Total number of genes in the category
        M : int
            Total number of genes
        """
        
        pval = hypergeom.sf(k - 1, M, n, N, loc=0)
        
        return pval
    
    ### ------------------------------------ ###
    
    @staticmethod
    def fisher_testing(aa, ab, ba, bb):
        
        """
        Hypergeometric test for enrichment
        
        Parameters
        ----------
        aa : int
            Number of SIGNIFICANT genes IN the category
        ab : int
            Number of SIGNIFICANT genes NOT IN the category
        ba : int
            Number of NOT SIGNIFICANT genes IN the category
        bb : int
            Number of NOT SIGNIFICANT genes NOT IN the category
        """
        
        table = [[aa, ab],
                 [ba, bb]]
        
        pval = fisher_exact(table, alternative='greater')[1]
        
        return pval
    
    ### ------------------------------------ ###
    ### VISUALIZATION                        ###
    ### ------------------------------------ ###
    
    def plot_data(self, p_thr=0.05, max_categories=15):
        
        """
        Class for the over-representation or depletion analysis of gene sets.

        Parameters
        ----------
        p_thr : float, optional
            Threshold for raw p-value enrichment/depletion.
            Determines how many categories to plot.
            Default = 0.05
        max_categories : int, optional
            Maximum number of categories to plot.
            Default = 15
        """
        
        # Filter results
        data = self.results.copy()
        data = data.loc[data.padj < p_thr,].iloc[:15,]

        # Change Inf and -Inf values
        vals = data.enrichment_score[data.enrichment_score.abs() != np.inf]
        data.replace(np.inf, max(1, vals.max()), inplace=True)
        data.replace(-np.inf, min(-1, vals.min()), inplace=True)
        
        if data.shape[0] >= 3: # Only plot if there's at least 3 categories
        
            # Make category names more readable
            data.category = [self.break_string(c) for c in data.category.values]
            
            # Add log10(- padj) and number of hits
            data['-log10(padj)'] = - np.log10(data.padj.values)
            data['Count'] = [int(d.gene_ratio.split('/')[0]) for _,d in data.iterrows()]
            
            # Sort data by enrichment_score
            data.sort_values(by='enrichment_score', inplace=True, ascending=False)
            
            # Set palette for -log10(padj)
            palette = sns.light_palette('red', n_colors=50, reverse=False, as_cmap=True)
            
            # Check if data has to be split into enriched/depleted
            if len(set(data.behaviour)) == 2:
                
                if (data.behaviour == 'enrichment').sum() == (data.behaviour == 'depletion').sum():
                    
                    width_ratios = [1, 1]
                
                elif (data.behaviour == 'enrichment').sum() > (data.behaviour == 'depletion').sum():
                    
                    width_ratios = [0.5, 1]
                
                else:
                    
                    width_ratios = [1, 0.5]
                
                figsize = (10, data.shape[0] // 2)
                
                fig, (ax1, ax2) = plt.subplots(nrows=1,
                                               ncols=2,
                                               figsize=figsize,
                                               sharex=False,
                                               sharey=True,
                                               width_ratios=width_ratios)
                
                hue_norm = (data['-log10(padj)'].min(), data['-log10(padj)'].max())

                data_depletion = data.copy()
                data_depletion.loc[data.behaviour == 'enrichment', 'enrichment_score'] = np.nan
                sns.scatterplot(ax=ax1, data=data_depletion, x='enrichment_score', y='category', hue='-log10(padj)', hue_norm=hue_norm, size='Count',  sizes=(10,100), palette=palette, edgecolor='black', lw=1)
                ax1.set_title('Depleted categories')
                ax1.set_xlim(data_depletion.enrichment_score.min() * 1.5, data_depletion.enrichment_score.max() * 0.5)
                ax1.get_legend().remove()
                ax1.set_xlabel('')
                ax1.set_ylabel('')
                
                data_enrichment = data.copy()
                data_enrichment.loc[data.behaviour == 'depletion', 'enrichment_score'] = np.nan
                sns.scatterplot(ax=ax2, data=data_enrichment, x='enrichment_score', y='category', hue='-log10(padj)', hue_norm=hue_norm, size='Count', sizes=(10, 100), palette=palette, edgecolor='black', lw=1)
                ax2.set_title('Enriched categories')
                ax2.set_xlim(data_enrichment.enrichment_score.min() * 0.5, data_enrichment.enrichment_score.max() * 1.5)
                ax2.set_xlabel('')
                ax2.set_ylabel('')
                
                sns.move_legend(ax2, "upper left", bbox_to_anchor=(1, 1))
                plt.subplots_adjust(wspace=0.05, hspace=0)
                
                fig.text(0, -0.05, 'Enrichment score', ha='center', va='top', transform=ax2.transAxes)
    
                plt.tight_layout()
                
                plt.savefig(f'{self.output_prefix}.png', dpi=300)
                plt.close()
            
            else:
                
                figsize = (10, data.shape[0] // 2)
                
                plt.figure(figsize=figsize)
                
                ax = sns.scatterplot(data, x='enrichment_score', y='category', hue='-log10(padj)', size='Count', sizes=(10, 100), palette=palette, lw=1, edgecolor='black')
                
                sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
                plt.ylabel('')
                plt.xlabel('Enrichment score')
                plt.ylabel('')
    
                plt.tight_layout()
                
                plt.savefig(f'{self.output_prefix}.png', dpi=300)
                plt.close()
    
    ### ------------------------------------ ###
    
    @staticmethod
    def break_string(string, max_chars_per_line=30):
        
        """
        Break the category names to be more readable.
        """
        
        words = string.split(' ')
        
        split_txt = ''
        breaks = 1
        for w in words:
          
          if len(split_txt + ' ' + w) / max_chars_per_line > breaks:
            
            breaks += 1
            split_txt = split_txt + '\n' + w
            
          elif len(split_txt) != 0:
            
            split_txt = split_txt + ' ' + w
            
          else:
            
            split_txt = w
        
        return split_txt
        
### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection
from sys import argv

if '--help' in argv:
    
    print_help()
    exit()

else:
    
    diff_analysis, term2gene, universe, min_genes, max_genes, filter_genes, output_prefix, split_updown, enrichment_type = parse_args()

if split_updown:
    
    # All genes
    gene_set = diff_analysis.gene_symbol.values.copy()
    analysis = enrichment_analysis(gene_set, term2gene, universe, f'{output_prefix}_AllGenes', filter_genes)
    analysis.process_data(min_genes=min_genes, max_genes=max_genes, enrichment_type=enrichment_type)
    analysis.plot_data(p_thr=0.05, max_categories=15)
    
    # Upreg only
    gene_set = diff_analysis.loc[diff_analysis.log2FoldChange.values > 0, 'gene_symbol'].values.copy()
    analysis = enrichment_analysis(gene_set, term2gene, universe, f'{output_prefix}_UpregGenes', filter_genes)
    analysis.process_data(min_genes=min_genes, max_genes=max_genes, enrichment_type=enrichment_type)
    analysis.plot_data(p_thr=0.05, max_categories=15)
    
    # Downreg only
    gene_set = diff_analysis.loc[diff_analysis.log2FoldChange.values < 0, 'gene_symbol'].values.copy()
    analysis = enrichment_analysis(gene_set, term2gene, universe, f'{output_prefix}_DownregGenes', filter_genes)
    analysis.process_data(min_genes=min_genes, max_genes=max_genes, enrichment_type=enrichment_type)
    analysis.plot_data(p_thr=0.05, max_categories=15)

else:
    
    gene_set = diff_analysis.gene_symbol.values.copy()
    analysis = enrichment_analysis(gene_set, term2gene, universe, output_prefix, filter_genes)
    analysis.process_data(min_genes=min_genes, max_genes=max_genes, enrichment_type=enrichment_type)
    analysis.plot_data(p_thr=0.05, max_categories=15)
