#!/usr/bin/env python3

"""
Interpolation of MEME SEA and ATAC Seq to find genes with TFBSs for a specific TF
"""

### ---------------------------------------- ###

def parse_args():

    print(argv)
    
    tf = argv[argv.index('--tf') + 1]
    tf = tf.upper()
    
    meme_sea_file = argv[argv.index('--sea') + 1]
    meme_sea = pd.read_csv(meme_sea_file, sep='\t')
    
    meme_sea_seqs_file = argv[argv.index('--sea_seqs') + 1]
    meme_sea_seqs = pd.read_csv(meme_sea_seqs_file, sep='\t')
    
    atac_file = argv[argv.index('--atac') + 1]
    atac = pd.read_csv(atac_file, sep='\t')
    
    return tf, meme_sea, meme_sea_seqs, atac

### ------------------MAIN------------------ ###

import pandas as pd

from sys import argv

# Parse args
tf, meme_sea, meme_sea_seqs, atac = parse_args()

if tf in meme_sea.ALT_ID.values:
    
    # Find optimal threshold value for sea score
    thr = meme_sea.loc[meme_sea.ALT_ID == tf, 'SCORE_THR'].values[0]
    
    # Filter meme_sea and extract peaks ids (holdhout seqs are discarded)
    peaks_ids = meme_sea_seqs.loc[(meme_sea_seqs.motif_ALT_ID == tf) &
                                  (meme_sea_seqs.seq_Score >= thr) &
                                  (meme_sea_seqs.seq_Class == 'tp') &
                                  (meme_sea_seqs['is_holdout?'] == 0),
                                  'seq_ID'].values
    peaks_ids = [int(pid.replace('peak_', '')) for pid in peaks_ids]
    
    # Extract list of genes that contributed to the enrichment of the TF
    genes_with_binding_site = atac.loc[atac.PeakID.isin(peaks_ids), ["GeneID", "GeneName"]]
    genes_with_binding_site.drop_duplicates(inplace=True)
    
    # Write to file
    genes_with_binding_site.to_csv(f'genes_bound_by_{tf}.tsv', sep='\t', index=False)
    