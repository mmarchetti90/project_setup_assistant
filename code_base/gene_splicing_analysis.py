#!/usr/bin/env python3

"""
Given a set of BAM files, target genes, and a gtf file, the class finds reads with junctions across exons
Data can be visualized for separate samples or as an average of selected samples
"""

### IMPORTS -------------------------------- ###

import numpy as np
import pandas as pd
import pysam
import re
import seaborn as sns

from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.stats import linregress, mannwhitneyu
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families.family import NegativeBinomial
from statsmodels.stats.multitest import fdrcorrection

### CLASSES AND FUNCTIONS ------------------ ###

class splice_check:

    def __init__(self, genes, gtf_path, bam_paths, sample_ids):

        # List of genes of interest
        
        self.genes = genes

        # Exon coordinates of genes of interest
        
        self.exon_coords = self.parse_gtf(gtf_path, genes)

        # Sample IDs
        
        self.sample_ids = sample_ids

        # Paths of bam files
        
        self.bam_files = {s : bp for s,bp in zip(sample_ids, bam_paths)}

    ### ------------------------------------ ###
    ### GTF PARSING                          ###
    ### ------------------------------------ ###

    @staticmethod
    def parse_gtf(path, target_genes=[]):
        
        """
        This function parses a GTF file and returns the exon coordinates of genes of interest
        N.B. If an exon has multiple IDs (e.g. belongs to different transcripts), then only the first ID is kept
        """
        
        print('### Parsing GTF file')
        
        # Load GTF
        
        gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
        gtf_data.columns = ['contig', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

        ### Only keep exons

        gtf_data = gtf_data.loc[gtf_data.feature == 'exon', ['contig', 'start', 'end', 'strand', 'attribute']]
        
        # Get biotype and gene id, gene symbol, transcript id, and exon id

        biotypes, gene_ids, gene_symbols, transcript_ids, exon_ids = [], [], [], [], []
        for _,row in gtf_data.iterrows():
            
            info = row.values[-1]
            
            biotype = re.findall('gene_biotype "\w+";', info)[0]
            biotype = biotype.replace('gene_biotype ', '').replace(';', '').replace('"', '')
            
            gene = re.findall('gene_id "\w+";', info)[0]
            gene = gene.replace('gene_id ', '').replace(';', '').replace('"', '')
            
            if 'gene_name' in info:
                
                symbol = info[info.index('gene_name "') + len('gene_name "'):]
                symbol = symbol[:symbol.index('"')]
            
            else:
                
                symbol = ''
            
            if 'transcript_id' in info:
                
                transcript = info[info.index('transcript_id "') + len('transcript_id "'):]
                transcript = transcript[:transcript.index('"')]
            
            else:
                
                transcript = ''
            
            if 'exon_id' in info:
                
                exon = info[info.index('exon_id "') + len('exon_id "'):]
                exon = exon[:exon.index('"')]
            
            else:
                
                exon = ''
            
            biotypes.append(biotype)
            gene_ids.append(gene)
            gene_symbols.append(symbol)
            transcript_ids.append(transcript)
            exon_ids.append(exon)

        gtf_data['biotype'] = biotypes
        gtf_data['gene_id'] = gene_ids
        gtf_data['gene_symbol'] = gene_symbols
        gtf_data['transcript_id'] = transcript_ids
        gtf_data['exon_id'] = exon_ids
        
        # Only keep genes of interest
        
        gtf_data = gtf_data.loc[(gtf_data['gene_id'].isin(target_genes)) |
                                (gtf_data['gene_symbol'].isin(target_genes)),]
        
        # Fix dtypes
        
        gtf_data.loc[:, ['start', 'end']] = gtf_data[['start', 'end']].astype(int)
        
        # Drop duplicate exons
        
        gtf_data = gtf_data.drop_duplicates(subset=['gene_id', 'contig', 'start', 'end'], keep='first')

        # Rename exons

        for g in target_genes:

            g_pos = np.where((gtf_data['gene_id'] == g) | (gtf_data['gene_symbol'] == g))[0]
            
            g_strand = gtf_data.iloc[g_pos[0]]['strand']
            
            if g_strand != '-':
                
                g_exons_starts = gtf_data.iloc[g_pos]['start'].values
                
                sort_order = np.argsort(g_exons_starts)
            
            else:
                
                g_exons_starts = gtf_data.iloc[g_pos]['end'].values
                
                sort_order = np.argsort(g_exons_starts)[::-1]
            
            new_exons_ids = ['' for _ in range(g_exons_starts.shape[0])]
            
            for n,o in enumerate(sort_order):
                
                new_exons_ids[o] = f'E{n+1}'
            
            exon_id_col_idx = gtf_data.columns.tolist().index('exon_id')
            
            gtf_data.iloc[g_pos, exon_id_col_idx] = new_exons_ids
        
        return gtf_data

    ### ------------------------------------ ###
    ### PARSE READS                          ###
    ### ------------------------------------ ###

    def parse_gene_reads(self, min_mapq=20, norm_target=1e6):

        """
        Wrapper for parse_sample_reads
        """

        # Init results dicts
        self.junctions_dict = {}
        self.coverages_dict = {}
        self.reads_stats = {}

        # Process genes
        
        for gene in self.genes:
            
            print(f'### Parsing reads for {gene}')

            self.junctions_dict[gene] = pd.DataFrame({col : []
                                                      for col in ['sample', 'start', 'stop', 'reads_count', 'norm_reads_count', 'left_exon', 'right_exon']})
            
            self.coverages_dict[gene] = pd.DataFrame({col : []
                                                      for col in ['sample', 'pos', 'cov']})

            for n,(sample,bam_path) in enumerate(self.bam_files.items()):
                
                print(f'Processing {sample}')

                junctions, coverage = self.parse_sample_reads(gene, sample, bam_path, min_mapq, norm_target)
                
                junctions.loc[:, 'sample'] = [sample for _ in range(junctions.shape[0])]
                junctions = junctions[['sample', 'start', 'stop', 'reads_count', 'norm_reads_count', 'left_exon', 'right_exon']]
                coverage = pd.DataFrame({'sample' : [sample for _ in range(coverage.shape[1])],
                                         'pos' : coverage[0,],
                                         'cov' : coverage[1,]})

                self.junctions_dict[gene] = pd.concat([self.junctions_dict[gene], junctions.copy()], axis=0, ignore_index=True)
                self.coverages_dict[gene] = pd.concat([self.coverages_dict[gene], coverage.copy()], axis=0, ignore_index=True)

            # Save data to file
            
            self.junctions_dict[gene].to_csv(f'{gene}_junctions.tsv', sep='\t', index=False)
            
            self.coverages_dict[gene].to_csv(f'{gene}_coverage.tsv', sep='\t', index=False)

    ### ------------------------------------ ###

    def parse_sample_reads(self, gene, sample_name, bam_path, min_mapq=20, norm_target=1e6):

        """
        Parsing aligned reads to find gene coverage and splicing events
        """

        # Get exons for specific gene
        
        gene_exon_coords = self.exon_coords.loc[(self.exon_coords.gene_id == gene) |
                                                (self.exon_coords.gene_symbol == gene),].copy()
        
        gene_contig = gene_exon_coords['contig'].values[0]
        gene_start = gene_exon_coords['start'].min()
        gene_end = gene_exon_coords['end'].max()
        
        # Get number of reads mapping to the gene, then compute a normalization factor for read counts
        
        tot_gene_reads = self.get_gene_counts(bam_path, gene_contig, gene_start, gene_end, min_mapq)
        
        norm_factor = norm_target / tot_gene_reads
        
        # Get introns, i.e. splice sites
        
        unique_junctions = self.get_junctions(bam_path, gene_exon_coords, min_mapq, norm_factor)
        
        # Compute coverage
            
        gene_cov = self.get_coverage(bam_path, gene_contig, gene_start, gene_end, min_mapq)

        return unique_junctions, gene_cov
    
    ### ------------------------------------ ###
    ### COMPARE JUNCTIONS DATA               ###
    ### ------------------------------------ ###
    
    def compare_gene_junctions(self, gene, samples_1, samples_2, group_1_name='g1', group_2_name='g2', collapse_replicates=False, canonical_only=False, kneedle_filter=False, min_reads=1):

        """
        Comparing gene splicing junctions between two samples/groups
        Different analyses depending on whether or not replicates are not available for both conditions
        """
        
        if collapse_replicates or len(samples_1) == 1 or len(samples_2) == 1:
            
            comparison_data = self.single_sample_comparison(gene, samples_1, samples_2, group_1_name, group_2_name, canonical_only, kneedle_filter, min_reads)
        
        else:
            
            comparison_data = self.multi_sample_comparison(gene, samples_1, samples_2, group_1_name, group_2_name, canonical_only, min_reads)
        
        # Write to file
        
        comparison_data.to_csv(f'{gene}_{group_1_name}_vs_{group_2_name}_junctions.tsv', sep='\t', index=False)
        
        return comparison_data
    
    ### ------------------------------------ ###
    
    def single_sample_comparison(self, gene, samples_1, samples_2, group_1_name='g1', group_2_name='g2', canonical_only=False, kneedle_filter=False, min_reads=1):
        
        # Get samples_1 junctions
        
        #sample_1_name = samples_1[0] if len(samples_1) == 1 else group_1_name
        sample_1_name = group_1_name
        
        u_junc_1 = self.collapse_group_junctions(samples_1, self.junctions_dict[gene].copy(), canonical_only)
        u_junc_1 = u_junc_1[['start', 'stop', 'left_exon', 'right_exon', 'reads_count', 'norm_reads_count']]
        u_junc_1.columns = ['start', 'stop', 'left_exon', 'right_exon', f'reads_count_{sample_1_name}', f'norm_reads_count_{sample_1_name}']
        
        # Get samples_2 junctions
        
        #sample_2_name = samples_2[0] if len(samples_2) == 1 else group_2_name
        sample_2_name = group_2_name
        
        u_junc_2 = self.collapse_group_junctions(samples_2, self.junctions_dict[gene].copy(), canonical_only)
        u_junc_2 = u_junc_2.loc[u_junc_2['reads_count'] > min_reads,]
        u_junc_2.columns = ['start', 'stop', 'left_exon', 'right_exon', f'reads_count_{sample_2_name}', f'norm_reads_count_{sample_2_name}']
        
        # Merge
        
        all_data = pd.merge(u_junc_1, u_junc_2,
                            how='outer', on=['start', 'stop', 'left_exon', 'right_exon'])
        
        # Fill missing count values

        all_data[all_data.columns[4:]] = all_data[all_data.columns[4:]].fillna(value=0)
        
        # Filter junctions based on raw reads
        
        row_filter = (all_data[[f'reads_count_{s}'
                                for s in [sample_1_name, sample_2_name]]].values < min_reads).sum(axis=1) < 2

        all_data = all_data.loc[row_filter,]
        
        # Filter junctions with Kneedle
        
        if kneedle_filter:
        
            thresholds = {s : self.kneedle(all_data[f'norm_reads_count_{s}'].values)[1] for s in [sample_1_name, sample_2_name]}
            row_filter = (all_data[[f'norm_reads_count_{s}'
                                    for s in [sample_1_name, sample_2_name]]].values < list(thresholds.values())).sum(axis=1) < 2
        
            all_data = all_data.loc[row_filter,]
        
        else:
            
            thresholds = {s : 0 for s in [sample_1_name, sample_2_name]}
        
        ### Reset index
        
        all_data = all_data.reset_index(drop=True)
        
        ### Tag junction
        
        # Find junctions unique to sample_1_name
        
        new_col_name = f'{sample_1_name}_unique'
        new_col_data = np.array(['' for _ in range(all_data.shape[0])], dtype='<U64')
        row_filter = ((all_data[f'norm_reads_count_{sample_1_name}'] > thresholds[sample_1_name]) &
                      (all_data[f'norm_reads_count_{sample_2_name}'] <= thresholds[sample_2_name]))
        new_col_data[row_filter] = 'FLAG'
        all_data.loc[:, new_col_name] = new_col_data
        
        # Find junctions missing in sample_2_name
        
        new_col_name = f'{sample_2_name}_unique'
        new_col_data = np.array(['' for _ in range(all_data.shape[0])], dtype='<U64')
        row_filter = ((all_data[f'norm_reads_count_{sample_2_name}'] > thresholds[sample_2_name]) &
                      (all_data[f'norm_reads_count_{sample_1_name}'] <= thresholds[sample_1_name]))
        new_col_data[row_filter] = 'FLAG'
        all_data.loc[:, new_col_name] = new_col_data
        
        # Find junctions that break correlation between samples
        
        new_col_name = f'{sample_1_name}_vs_{sample_2_name}'
        new_col_data = np.array(['' for _ in range(all_data.shape[0])], dtype='<U64')

        junctions_subset = all_data.loc[(all_data[f'norm_reads_count_{sample_1_name}'] > 0) & (all_data[f'norm_reads_count_{sample_2_name}'] > 0),]
        
        x_val = junctions_subset[f'norm_reads_count_{sample_1_name}'].values
        y_val = junctions_subset[f'norm_reads_count_{sample_2_name}'].values
        
        outliers_idx, r, p = self.find_outliers(x_val, y_val, sample_1_name, sample_2_name, f'{gene}_{sample_1_name}_vs_{sample_2_name}_linregress.png')
        outliers_idx = junctions_subset.index[outliers_idx]
        
        new_col_data[outliers_idx] = 'FLAG'
        all_data.loc[:, new_col_name] = new_col_data
        
        # Plot flagged junctions
        
        if 'FLAG' in all_data.iloc[:, -3:].values:
            
            self.plot_flagged_junctions(gene, all_data, f'{gene}_{sample_1_name}_vs_{sample_2_name}_junctions.png')
        
        return all_data
    
    ### ------------------------------------ ###
    
    def multi_sample_comparison(self, gene, samples_1, samples_2, group_1_name='g1', group_2_name='g2', canonical_only=False, min_reads=1):
        
        # Subset junctions
        
        junc = self.junctions_dict[gene].copy()
        
        junc = junc.loc[junc['sample'].isin(samples_1 + samples_2),]
        
        if canonical_only:
            
            junc = junc.loc[(junc['left_exon'] != '') |
                            (junc['right_exon'] != ''),]
        
        # Statistics
        
        N = 8
        
        stats = []
        stats_header = ['start', 'stop',
                        f'{group_1_name}_norm_median', f'{group_2_name}_norm_median',
                        'stat_type', 'stat', 'pval', 'padj']
        
        for start,stop in junc[['start', 'stop']].drop_duplicates(keep='first').values:
            
            s1_dat = junc.loc[(junc['start'] == start) &
                              (junc['stop'] == stop) &
                              (junc['sample'].isin(samples_1)),
                              'norm_reads_count'].values

            s1_dat = np.concatenate([s1_dat, np.zeros(max(0, len(samples_1) - s1_dat.shape[0]))]) # Add missing samples as undetected
        
            s2_dat = junc.loc[(junc['start'] == start) &
                              (junc['stop'] == stop) &
                              (junc['sample'].isin(samples_2)),
                              'norm_reads_count'].values

            s2_dat = np.concatenate([s2_dat, np.zeros(max(0, len(samples_2) - s2_dat.shape[0]))]) # Add missing samples as undetected

            s1_median = np.median(s1_dat)
            
            s2_median = np.median(s2_dat)
            
            if (s1_dat > 0).sum() < 3 and (s2_dat > 0).sum() < 3:
                
                continue
            
            elif (s1_dat > 0).sum() < 3:
                
                median_reads = junc.loc[(junc['start'] == start) &
                                        (junc['stop'] == stop) &
                                        (junc['sample'].isin(samples_2)),
                                        'reads_count'].values

                median_reads = np.concatenate([median_reads, np.zeros(max(0, len(samples_2) - median_reads.shape[0]))])

                median_reads = np.median(median_reads)
                
                if median_reads >= min_reads:
                
                    stats.append([start, stop, s1_median, s2_median, f'{group_2_name}_unique', 0, 1, 0])
                    
                else:
                    
                    continue
            
            elif (s2_dat > 0).sum() < 3:
                
                median_reads = junc.loc[(junc['start'] == start) &
                                        (junc['stop'] == stop) &
                                        (junc['sample'].isin(samples_1)),
                                        'reads_count'].values

                median_reads = np.concatenate([median_reads, np.zeros(max(0, len(samples_1) - median_reads.shape[0]))])

                median_reads = np.median(median_reads)
                
                if median_reads >= min_reads:
                
                    stats.append([start, stop, s1_median, s2_median, f'{group_1_name}_unique', 0, 1, 0])
                    
                else:
                    
                    continue
                
            elif s1_dat.shape[0] >= N and s2_dat.shape[0] >= N:
                
                # Mann-Whitney U test
                
                stat_type = 'utest'
                
                stat_val = np.log2(s1_median / (s2_median + 1e-12))
                
                pval = mannwhitneyu(s1_dat, s2_dat, alternative='two-sided').pvalue
            
            else:
                
                # Negative binomial + Wald test
                
                stat_type = 'glm_nb'
                
                glm_data = pd.DataFrame({'junction_counts' : np.concatenate([s1_dat, s2_dat]),
                                         'sample_group' : [group_1_name for _ in range(s1_dat.shape[0])] + [group_2_name for _ in range(s2_dat.shape[0])]})
                
                model = GLM.from_formula('junction_counts ~ C(sample_group)', glm_data, family=NegativeBinomial())
                
                res = model.fit()
                
                stat_val = [st for param,st in res.params.items() if 'sample_group' in param][0]
                
                pval = [p for param,p in res.pvalues.items() if 'sample_group' in param][0]
        
            stats.append([start, stop, s1_median, s2_median, stat_type, stat_val, pval, 0.])
        
        stats = pd.DataFrame(stats, columns=stats_header)
        
        stats.loc[:, 'padj'] = fdrcorrection(stats['pval'].values, alpha=0.05, is_sorted=False)[1]
        
        return stats
    
    ### ------------------------------------ ###
    
    @staticmethod
    def find_outliers(x_val, y_val, x_name, y_name, diagnostic_plot_name='linregress.png'):

        """
        Assuming that the proband and a relative use mostly the same splicing junctions in a
        similar way, then the normalized counts for shared junctions should correlate
        A linear regression is fit and outliers are defined as being > 2 standard deviations from
        the fit
        """
        
        # Stats
        slope, intercept, r, p, std_err = linregress(x_val, y_val, alternative='two-sided')
        std = std_err * (len(x_val))**0.5
        y_predict = np.array([slope * x + intercept for x in x_val])
        
        # Define outliers based on > 2 std from best-fit line
        outliers_idx = np.where(abs(y_val - y_predict) > 2 * std)[0]
        
        if len(diagnostic_plot_name):
        
            # Plot data
            plt.figure(figsize=(5, 5))
            plt.title(f'r = {r:.3f}\np = {p:.2e}', loc='left')
            plt.xlabel(f'{x_name} counts')
            plt.ylabel(f'{y_name} counts')
            plt.scatter(x_val, y_val, color='lightgray', lw=0.5, edgecolor='black')
            plt.scatter(x_val[outliers_idx], y_val[outliers_idx], color='red', lw=0.5, edgecolor='black')
            plt.plot([x_val.min(), x_val.max()], [y_predict.min(), y_predict.max()], 'green')
            plt.plot([x_val.min(), x_val.max()], [y_predict.min() - 2 * std, y_predict.max() - 2 * std], 'black', linestyle='--')
            plt.plot([x_val.min(), x_val.max()], [y_predict.min() + 2 * std, y_predict.max() + 2 * std], 'black', linestyle='--')
            plt.tight_layout()
            plt.savefig(diagnostic_plot_name, dpi=300)
            plt.close()
        
        return outliers_idx, r, p

    ### ------------------------------------ ###
    ### UTILS                                ###
    ### ------------------------------------ ###
    
    @staticmethod
    def collapse_group_coverages(group_samples, cov):
        
        cov = cov.loc[cov['sample'].isin(group_samples),]
        
        if len(group_samples) == 1:
            
            cov = cov.drop(columns='sample')
        
        else:
            
            cov = cov.groupby(by=['pos'])['cov'].sum().reset_index(drop=False)
            
            cov.loc[:, 'cov'] *= (1 / len(group_samples))
            
            cov = cov.sort_values('pos')
        
        return cov
        
    ### ------------------------------------ ###
    
    @staticmethod
    def collapse_group_junctions(group_samples, junc, canonical_only=True):
        
        junc = junc.loc[junc['sample'].isin(group_samples),]
        
        if len(group_samples) == 1:
            
            junc = junc.drop(columns='sample')
        
        else:
            
            junc = junc.groupby(by=['start', 'stop', 'left_exon', 'right_exon'])[['reads_count', 'norm_reads_count']].median().reset_index(drop=False)
        
        # Only keep junctions with known exons
        
        if canonical_only:
            
            #junc = junc.loc[(junc['left_exon'] != '') |
            #                (junc['right_exon'] != ''),]

            junc = junc.loc[(junc['left_exon'] != '') &
                            (junc['right_exon'] != ''),]
        
        return junc
    
    ### ------------------------------------ ###
    
    @staticmethod
    def get_coverage(bam_path, gene_contig, gene_start, gene_end, min_mapq=20):
        
        alignments = pysam.AlignmentFile(bam_path, 'rb')
        
        # Compute coverage
        
        try:
        
            gene_cov = alignments.count_coverage(contig=gene_contig, start=gene_start - 1, stop=gene_end, quality_threshold=min_mapq)
            gene_cov = np.max(gene_cov, axis=0)
        
        except:
            
            gene_cov = alignments.count_coverage(contig='chr'+gene_contig, start=gene_start - 1, stop=gene_end, quality_threshold=min_mapq)
            gene_cov = np.max(gene_cov, axis=0)
            
        gene_cov = np.stack([np.arange(gene_start, gene_end + 1, 1), gene_cov])
        
        return gene_cov
    
    ### ------------------------------------ ###
    
    @staticmethod
    def get_gene_counts(bam_path, gene_contig, gene_start, gene_end, min_mapq=20):
        
        # Open bam file
        
        alignments = pysam.AlignmentFile(bam_path, 'rb')
        
        try:
            
            reads_of_interest = alignments.fetch(contig=gene_contig, start=gene_start - 1, stop=gene_end)
        
        except:
            
            reads_of_interest = alignments.fetch(contig='chr' + gene_contig, start=gene_start - 1, stop=gene_end)
        
        # Filter by mapq
        
        filtered_reads_count = sum(1 for read in reads_of_interest if read.mapping_quality >= min_mapq and read.is_mapped)
        
        return filtered_reads_count
    
    ### ------------------------------------ ###
    
    @staticmethod
    def get_junctions(bam_path, ex_coords, min_mapq=20, norm_factor=1):
        
        # Open bam file
        
        alignments = pysam.AlignmentFile(bam_path, 'rb')

        # Fetch reads
        
        gene_contig = ex_coords['contig'].values[0]
        gene_start = ex_coords['start'].min()
        gene_end = ex_coords['end'].max()
        
        try:
            
            reads_of_interest = alignments.fetch(contig=gene_contig, start=gene_start - 1, stop=gene_end)
        
        except:
            
            reads_of_interest = alignments.fetch(contig='chr' + gene_contig, start=gene_start - 1, stop=gene_end)
        
        # Filter by mapq
        
        filtered_reads = (read for read in reads_of_interest if read.mapping_quality >= min_mapq and read.is_proper_pair)
        
        # Get introns, i.e. splice sites
        
        introns = alignments.find_introns(filtered_reads)
        
        unique_junctions = []
        unique_junctions_header = ['start', 'stop', 'reads_count', 'norm_reads_count', 'left_exon', 'right_exon']
        
        for (start, stop), count in introns.items():
            
            stop += 1
            
            norm_count = count * norm_factor
            
            left = ','.join(ex_coords.loc[ex_coords['end'] == start, 'exon_id'].values.astype(str))
            right = ','.join(ex_coords.loc[ex_coords['start'] == stop, 'exon_id'].values.astype(str))
            
            unique_junctions.append([start, stop, count, norm_count, left, right])
        
        unique_junctions = pd.DataFrame(unique_junctions, columns=unique_junctions_header)
        
        return unique_junctions
    
    ### ------------------------------------ ###
    
    @staticmethod
    def kneedle(vector):
    
        """
        Kneedle to find threshold cutoff.
        """
        
        # Sort data
        vector = np.sort(vector)[::-1]
        
        # Find gradient and intercept
        x0, x1 = 0, len(vector)
        y0, y1 = max(vector), min(vector)
        gradient = (y1 - y0) / (x1 - x0)
        intercept = y0
        
        # Compute difference vector
        difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(vector)]
        
        # Find max of difference_vector and define cutoff
        cutoff_index = difference_vector.index(max(difference_vector))
        cutoff_value = vector[cutoff_index]
        
        return cutoff_index, cutoff_value
    
    ### ------------------------------------ ###
    ### PLOTTING                             ###
    ### ------------------------------------ ###
    
    def plot_paired_sashimi(self, gene, samples_1, samples_2, group_1_name='', group_2_name='', canonical_only=False, min_reads=1, out_name='sashimi.png', color_1='red', color_2='deepskyblue'):
        
        colors = {}
        
        # Get samples_1 coverage and junctions
        
        sample_1_name = samples_1[0] if len(samples_1) == 1 else group_1_name
        
        colors[sample_1_name] = color_1
        
        g_cov_1 = self.collapse_group_coverages(samples_1, self.coverages_dict[gene].copy())
        
        u_junc_1 = self.collapse_group_junctions(samples_1, self.junctions_dict[gene].copy(), canonical_only)
        u_junc_1 = u_junc_1.loc[u_junc_1['reads_count'] >= min_reads,]
        u_junc_1 = u_junc_1[['start', 'stop', 'left_exon', 'right_exon', 'reads_count', 'norm_reads_count']]
        u_junc_1.columns = ['start', 'stop', 'left_exon', 'right_exon', f'reads_count_{sample_1_name}', f'norm_reads_count_{sample_1_name}']
        
        # Get samples_2 coverage and junctions
        
        sample_2_name = samples_2[0] if len(samples_2) == 1 else group_2_name
        
        colors[sample_2_name] = color_2
        
        g_cov_2 = self.collapse_group_coverages(samples_2, self.coverages_dict[gene].copy())
        
        u_junc_2 = self.collapse_group_junctions(samples_2, self.junctions_dict[gene].copy(), canonical_only)
        u_junc_2 = u_junc_2.loc[u_junc_2['reads_count'] >= min_reads,]
        u_junc_2 = u_junc_2[['start', 'stop', 'left_exon', 'right_exon', 'reads_count', 'norm_reads_count']]
        u_junc_2.columns = ['start', 'stop', 'left_exon', 'right_exon', f'reads_count_{sample_2_name}', f'norm_reads_count_{sample_2_name}']
        
        # Get exons coordinates
        
        ex_coords = self.exon_coords.loc[(self.exon_coords.gene_id == gene) |
                                         (self.exon_coords.gene_symbol == gene),].copy()
    
        # Aestetic parameters
        
        intron_thickness = 1
        exon_thickness = 8
        splicing_thickness = 1
        min_height = 1
        max_height = 10
        
        gene_start, gene_stop = ex_coords['start'].min(), ex_coords['end'].max()

        fig_width = max(10, abs(gene_stop - gene_start) / 10000)
        fig_height = 12
        
        # Add splicing size to u_junc, normalized
        
        juncs_plot_data = pd.DataFrame([[sample, start, stop, abs(stop - start)]
                                        for sample,j in zip([sample_1_name, sample_2_name], [u_junc_1, u_junc_2])
                                        for _,(start, stop, *_) in j.iterrows()],
                                       columns=['sample', 'start', 'stop', 'splicing_size'])
        
        juncs_plot_data = juncs_plot_data.assign(splicing_size_norm = min_height + np.log2(juncs_plot_data.splicing_size / juncs_plot_data.splicing_size.min()))
        
        # Normalize coverage to 0-1 range, then set g_cov_norm max height to 2 * juncs_plot_data.splicing_size_norm.max()

        g_cov_max = max([c['cov'].max() for c in [g_cov_1, g_cov_2]])
        
        g_cov_norm = {}
        
        for sample,cov in zip([sample_1_name, sample_2_name], [g_cov_1, g_cov_2]):
            
            g_cov_norm[sample] = cov.copy()
            
            g_cov_norm[sample].loc[:, 'cov'] = (g_cov_norm[sample]['cov'] / g_cov_max)
            g_cov_norm[sample].loc[:, 'cov'] = g_cov_norm[sample]['cov'] * (2 * juncs_plot_data.splicing_size_norm.max()) if len(juncs_plot_data) else max_height
        
        # Init plot
        
        plt.figure(figsize=(fig_width, fig_height))
        
        ### Gene structure
        
        # Plot gene length
        
        gene_structure_pos = -3
        plt.hlines(gene_structure_pos, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot exons
        
        for _, exon in ex_coords.iterrows():
            
            plt.hlines(gene_structure_pos, exon['start'], exon['end'], linestyles='solid', colors='black', linewidth=exon_thickness)
        
        ### Coverage and splicing
        
        offset = 0
        
        ordered_samples = list(g_cov_norm.keys())[::-1]
        
        for n,(s) in enumerate(ordered_samples):
            
            c = g_cov_norm[s].copy()
            
            j = juncs_plot_data.loc[juncs_plot_data['sample'] == s,]
            
            offset += 3
            
            offset += j.loc[((j.start >= gene_start) & (j['start'] <= gene_stop)) |
                            ((j.stop >= gene_start) & (j['stop'] <= gene_stop)), 'splicing_size_norm'].max()
        
            plt.hlines(offset, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
            
            # Plot coverage
            
            where_pos = [True] + list(np.diff(c['pos'].values) == 1) # Where to fill
            
            plt.fill_between(x=c['pos'].values, y1=np.zeros(c.shape[0]) + offset, y2=c['cov'].values + offset, where=where_pos, alpha=1, color=colors[s], linewidth=0)
            
            # Plot splicing events
            for _,splicing in j.iterrows():
                
                intercept = - splicing.splicing_size_norm
                
                if splicing.start < splicing.stop:
                    
                    spline = CubicSpline(x=[splicing['start'], (splicing['stop'] + splicing['start']) / 2, splicing['stop']],
                                         y=[0, intercept, 0])
                
                else:
                    
                    spline = CubicSpline(x=[splicing.stop, (splicing['stop'] + splicing['start']) / 2, splicing['start']],
                                         y=[0, intercept, 0])
                
                x = np.linspace(splicing['start'], splicing['stop'], 21)
                y = spline(x) + offset
                
                plt.plot(x, y, ls='-', lw=splicing_thickness, color='black', alpha=1)
            
            offset += c['cov'].max()
        
        plt.xlabel(f'{ex_coords["contig"].values[0]} (bp)')
        plt.xlim(gene_start, gene_stop)
        plt.yticks([])
        plt.title(None)
        plt.tight_layout()
        plt.box(False)
        
        plt.savefig(out_name, dpi=300)
        plt.close()

    ### ------------------------------------ ###
    
    def plot_flagged_junctions(self, gene, juncs, out_name='flagged_junctions.png'):
        
        """
        Sashimi plot for the flagged junctions of a specific gene
        """
        
        # Aestetic parameters
        
        intron_thickness = 1
        exon_thickness = 5
        splicing_thickness = 1
        min_height = 1
        
        # Get exons coordinates
        
        ex_coords = self.exon_coords.loc[(self.exon_coords.gene_id == gene) |
                                         (self.exon_coords.gene_symbol == gene),].copy()
    
        gene_start, gene_end = ex_coords['start'].min(), ex_coords['end'].max()
        
        fig_width = max(10, abs(gene_end - gene_start) / 10000)
        fig_height = 4
        
        # Add splicing size to juncs, normalized
        
        juncs_plot_data = pd.DataFrame([[start, stop, abs(stop - start)] for _,(start, stop, count, *_) in juncs.iterrows()], columns=['start', 'stop', 'splicing_size'])
        juncs_plot_data = juncs_plot_data.assign(splicing_size_norm = min_height + np.log2(juncs_plot_data.splicing_size / juncs_plot_data.splicing_size.min()))
        
        # Init plot
        
        plt.figure(figsize=(fig_width, fig_height))
        
        ### Gene structure
        
        gene_structure_pos = - 1
        
        # Plot gene length
        
        plt.hlines(gene_structure_pos, gene_start, gene_end, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot exons
        
        for _, exon in ex_coords.iterrows():
            
            plt.hlines(gene_structure_pos, exon['start'], exon['end'], linestyles='solid', colors='black', linewidth=exon_thickness)
        
        ### Coverage and splicing
        
        plt.hlines(0, gene_start, gene_end, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot splicing events
        
        for _,splicing in juncs_plot_data.iterrows():
            
            intercept = splicing.splicing_size_norm
            
            if splicing.start < splicing.stop:
                
                spline = CubicSpline(x=[splicing['start'], (splicing['stop'] + splicing['start']) / 2, splicing['stop']],
                                     y=[0, intercept, 0])
            
            else:
                
                spline = CubicSpline(x=[splicing['stop'], (splicing['stop'] + splicing['start']) / 2, splicing['start']],
                                     y=[0, intercept, 0])
            
            x = np.linspace(splicing['start'], splicing['stop'], 21)
            y = spline(x)
            
            plt.plot(x, y, ls='-', lw=splicing_thickness, color='black', alpha=1)
        
        ### FLAGS
        
        # Find columns with flags
        
        flags_cols = juncs.columns.values[-3:]
        
        # Find x and y axis positions as well as flag type
        
        base_y = gene_structure_pos = juncs_plot_data['splicing_size_norm'].max() + 2
        xvals, yvals, flags = [], [], []
        for _,j in juncs.iterrows():
            
            x = (j['start'] + j['stop']) / 2
            for n,col in enumerate(flags_cols):
                
                if j[col] == 'FLAG':
                    
                    y = base_y + n
                    
                    xvals.append(x)
                    yvals.append(y)
                    flags.append(col)
        
        # Convert to data frame
        
        flags_plot_data = pd.DataFrame(np.stack([xvals, yvals, flags], axis=-1), columns=['x', 'y', 'flags'])
        flags_plot_data = flags_plot_data.astype({'x' : float, 'y' : float, 'flags' : str})
        
        # Plotting flags
        
        scatter = sns.scatterplot(flags_plot_data, x='x', y='y', hue='flags', s=25, edgecolor='black', linewidth=1)
        sns.move_legend(scatter, loc='upper left', bbox_to_anchor=(1, 1), title='Flags')
        
        plt.xlabel(f'{ex_coords["contig"].values[0]} (bp)')
        plt.xlim(gene_start, gene_end)
        plt.ylabel(None)
        plt.yticks([])
        plt.title(None)
        plt.tight_layout()
        plt.box(False)
        
        plt.savefig(out_name, dpi=300)
        plt.close()
    
    ### ------------------------------------ ###

    def plot_single_sashimi(self, samples, gene, interval=[], canonical_only=True, min_reads=1, sample_group_name='', color='red'):
        
        """
        Sashimi and coverage plot for a specific gene and individual sample
        """
        
        # Get sample coverage and junctions
        
        g_cov = self.collapse_group_coverages(samples, self.coverages_dict[gene].copy())
        
        u_junc = self.collapse_group_junctions(samples, self.junctions_dict[gene].copy(), canonical_only)

        u_junc = u_junc.loc[u_junc['reads_count'] >= min_reads,]
        
        sample_name = samples[0] if len(samples) == 1 else sample_group_name
        
        # Get exons coordinates
        
        ex_coords = self.exon_coords.loc[(self.exon_coords.gene_id == gene) |
                                         (self.exon_coords.gene_symbol == gene),].copy()
    
        # Aestetic parameters
        
        intron_thickness = 1
        exon_thickness = 5
        splicing_thickness = 1
        min_height = 1
        max_height = 10
        
        if not len(interval):
            
            gene_start, gene_stop = ex_coords['start'].min(), ex_coords['end'].max()
            
        else:
            
            gene_start, gene_stop = interval
        
        fig_width = max(10, abs(gene_stop - gene_start) / 10000)
        fig_height = 4
        
        # Add splicing size to u_junc, normalized

        juncs_plot_data = pd.DataFrame([[start, end, abs(end - start)]
                                        for _,(start, end, *_) in u_junc.iterrows()],
                                       columns=['start', 'end', 'splicing_size'])

        juncs_plot_data.loc[:, 'splicing_size_norm'] = min_height + np.log2(juncs_plot_data['splicing_size'] / juncs_plot_data['splicing_size'].min())
        
        # Normalize coverage to 0-1 range
        
        g_cov.loc[:, 'cov'] = (g_cov['cov'] / g_cov['cov'].max())
        
        # Set g_cov max height to 2 * juncs_plot_data.splicing_size_norm.max()
        
        g_cov.loc[:, 'cov'] = g_cov['cov'] * (2 * juncs_plot_data['splicing_size_norm'].max()) if len(juncs_plot_data) else max_height
        
        # Init plot

        plt.figure(figsize=(fig_width, fig_height))
        
        ### Gene structure
        
        if len(juncs_plot_data):
            
            gene_structure_pos = - juncs_plot_data.loc[((juncs_plot_data['start'] >= gene_start) & (juncs_plot_data['start'] <= gene_stop)) |
                                                       ((juncs_plot_data['end'] >= gene_start) & (juncs_plot_data['end'] <= gene_stop)), 'splicing_size_norm'].max() - 2
        
        else:
            
            gene_structure_pos = - 1
        
        # Plot gene length
        
        plt.hlines(gene_structure_pos, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot exons
        
        for _, exon in ex_coords.iterrows():
            
            plt.hlines(gene_structure_pos, exon['start'], exon['end'], linestyles='solid', colors='black', linewidth=exon_thickness)
        
        ### Coverage and splicing
        
        plt.hlines(0, gene_start, gene_stop, linestyles='solid', colors='black', linewidth=intron_thickness)
        
        # Plot coverage
        
        where_pos = [True] + list(np.diff(g_cov['pos'].values) == 1) # Where to fill
        
        plt.fill_between(x=g_cov['pos'].values, y1=np.zeros(g_cov.shape[0]), y2=g_cov['cov'].values, where=where_pos, alpha=1, color=color, linewidth=0)
        
        # Plot splicing events
        
        for _,splicing in juncs_plot_data.iterrows():
            
            intercept = - splicing.splicing_size_norm
            
            if splicing['start'] < splicing['end']:
                
                spline = CubicSpline(x=[splicing['start'], (splicing['end'] + splicing['start']) / 2, splicing['end']],
                                     y=[0, intercept, 0])
            
            else:
                
                spline = CubicSpline(x=[splicing['end'], (splicing['end'] + splicing['start']) / 2, splicing['start']],
                                     y=[0, intercept, 0])
            
            x = np.linspace(splicing['start'], splicing['end'], 21)
            y = spline(x)
            
            plt.plot(x, y, ls='-', lw=splicing_thickness, color='black', alpha=1)
        
        plt.xlabel(f'{ex_coords.contig.values[0]} (bp)')
        plt.xlim(gene_start, gene_stop)
        plt.yticks([])
        plt.title(None)
        plt.tight_layout()
        plt.box(False)
        
        plt.savefig(f'{sample_name}_sashimi.png', dpi=300)
        plt.close()

### ---------------------------------------- ###

if __name__ == '__main__':
    
    pass
