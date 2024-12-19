#!/usr/bin/env python3

### ---------------------------------------- ###

def parse_args():

    # Experiment fimo file
    exp_fimo_file = argv[argv.index('--exp') + 1]

    # Number of sequences used for experiment
    exp_seqs = int(argv[argv.index('--exp_seqs') + 1])

    # Background fimo file
    bg_fimo_file = argv[argv.index('--bg') + 1]

    # Number of sequences used for background (assumed to contain also exp_seqs for the purpose of fisher_exact)
    bg_seqs = int(argv[argv.index('--bg_seqs') + 1])

    # p threshold
    if '--p_thr' in argv:

        p_thr = float(argv[argv.index('--p_thr') + 1])

    else:

        p_thr = 1e-4

    # Output prefix
    if '--out_name' in argv:

        out_name = argv[argv.index('--out_name') + 1]

    else:

        out_name = 'peaks_tf-summary'

    return exp_fimo_file, exp_seqs, bg_fimo_file, bg_seqs, p_thr, out_name

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from scipy.stats import fisher_exact
from sys import argv

### Parse args

exp_fimo_file, exp_seqs, bg_fimo_file, bg_seqs, p_thr, out_name = parse_args()

### Load data

# Load as Pandas dataframes
exp_fimo_data = pd.read_csv(exp_fimo_file, sep='\t', comment='#', header=0)
bg_fimo_data = pd.read_csv(bg_fimo_file, sep='\t', comment='#', header=0)

# For each gene, remove motifs at same coords
fields_subset = ['motif_alt_id', 'sequence_name', 'start', 'stop', 'strand']
exp_fimo_data = exp_fimo_file.drop_duplicates(subset=fields_subset, keep='first')
bg_fimo_data = bg_fimo_data.drop_duplicates(subset=fields_subset, keep='first')

# Filter for raw p-value
exp_fimo_data = exp_fimo_file.loc[exp_fimo_file.loc[:, 'p-value'].values < p_thr,]
bg_fimo_data = bg_fimo_data.loc[bg_fimo_data.loc[:, 'p-value'].values < p_thr,]

### Output the percentage of genes with possible binding sites with pval < p_thr

# Init summary table
results = ['\t'.join(['transcription_factor',
                      'experiment_sequences',
                      'background_sequences',
                      'experiment_sequences_with_possible_binding_sites',
                      'background_sequences_with_possible_binding_sites',
                      'frequency_in_experiment_sequences',
                      'frequency_in_background_sequences',
                      'fisher_test'])]

for tf in set(bg_fimo_data.motif_alt_id):
    
    # Subset experimental data
    exp_data_sub = exp_fimo_data.loc[exp_fimo_data.motif_alt_id.isin([tf])].copy()
    
    # Subset background data
    bg_data_sub = bg_fimo_data.loc[bg_fimo_data.motif_alt_id.isin([tf])].copy()

    # Get unique seqs
    exp_seqs_with_binding_sites = np.unique(exp_data_sub.sequence_name.values).shape[0]
    bg_seqs_with_binding_sites = np.unique(bg_data_sub.sequence_name.values).shape[0]
    
    # Get frequency
    exp_frequency = round(100 * exp_seqs_with_binding_sites / exp_seqs, 2)
    bg_frequency = round(100 * bg_seqs_with_binding_sites / bg_seqs, 2)
    
    # Fisher's test
    pval = fisher_exact([[exp_seqs_with_binding_sites, exp_seqs - exp_seqs_with_binding_sites],
                         [bg_seqs_with_binding_sites, bg_seqs - bg_seqs_with_binding_sites - exp_seqs]],
                        alternative='two-sided').pvalue

    # Add results
    results.append('\t'.join([tf,
                              str(exp_seqs),
                              str(bg_seqs),
                              str(exp_seqs_with_binding_sites),
                              str(bg_seqs_with_binding_sites),
                              str(exp_frequency),
                              str(bg_frequency),
                              str(pval)]))

results = '\n'.join(results)

with open(f'{out_name}.tsv', 'w') as output:
    
    output.write(results)