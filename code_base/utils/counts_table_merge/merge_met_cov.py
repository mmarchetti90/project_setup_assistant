#!/usr/bin/env python3

"""
Merges methylation coverage files and summarizes to cpg island level (if desired).
"""

### ---------------------------------------- ###

def parse_args():
    
    print(argv)
    
    # Coverage data folder
    cov_dir = argv[argv.index('--cov_dir') + 1]
    
    # Coverage file pattern
    pattern = argv[argv.index('--pattern') + 1]
    
    # Find files
    cov_files = {file.replace(pattern, '') : f'{cov_dir}/{file}' for file in listdir(cov_dir) if pattern in file}
    
    # CpG island coordinates
    if '--cpg_bed' in argv:
        
        cpg_bed_file = argv[argv.index('--cpg_bed') + 1]
        cpg_bed = pd.read_csv(cpg_bed_file,
                              sep='\t',
                              header=None,
                              names=['chrom', 'start', 'end', 'island_id'],
                              dtype={'chrom' : str,
                                     'start' : int,
                                     'end' : int,
                                     'island_id' : str})
        
        cpg_bed['chrom'] = [chrom.replace('chr', '') for chrom in cpg_bed['chrom'].values]
        
        cpg_bed.sort_values(by=['chrom', 'start'], inplace=True)
    
    else:
        
        cpg_bed = pd.DataFrame([], columns = ['chrom', 'start', 'end', 'island_id'])
    
    # Min depth (methylation of any CpG or island with depth below this threshold is set to 0)
    if'--min_depth' in argv:
        
        min_depth = int(argv[argv.index('--min_depth') + 1])
    
    else:
        
        min_depth = 50
        
    return cov_files, cpg_bed, min_depth

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from os import listdir
from sys import argv

### Parse args

cov_files, cpg_bed, min_depth = parse_args()

### Read and merge files

print('Merging files')

try:
    
    del all_data
        
except:
    
    pass

for file_id, file_path in cov_files.items():
    
    data = pd.read_csv(file_path, sep='\t',
                       header=None,
                       names=['chrom', 'start', 'end', 'methylation_percentage', 'count_methylated', 'count_unmethylated'],
                       dtype={'chrom' : str,
                              'start' : int,
                              'end' : int,
                              'methylation_percentage' : float,
                              'count_methylated' : int,
                              'count_unmethylated' : int})
    
    data.drop('methylation_percentage', axis=1, inplace=True)
    
    try:
        
        all_data = all_data.merge(data, on=['chrom', 'start', 'end'], how='outer')
        
        all_data.columns = np.concatenate([all_data.columns.values[:-2], data.columns.values[-2:] + f'_{file_id}'])
    
    except:
        
        all_data = data.copy()
        
        all_data.columns = ['chrom', 'start', 'end', f'count_methylated_{file_id}', f'count_unmethylated_{file_id}']

all_data.fillna(0, inplace=True)

all_data.sort_values(by=['chrom', 'start'], inplace=True)

### Merge to cytosine level (if desired)

if cpg_bed.shape[0]:
    
    print('Summarizing to CpG island level')
    
    methylation_data = []
    
    # Remove contigs not in common
    good_contigs = np.unique(cpg_bed.chrom.values)[np.isin(np.unique(cpg_bed.chrom.values), np.unique(all_data.chrom.values), assume_unique=True)]

    all_data = all_data.loc[all_data.chrom.isin(good_contigs),]
    cpg_bed = cpg_bed.loc[cpg_bed.chrom.isin(good_contigs),]
    
    for chrom in good_contigs:
        
        all_data_sub = all_data.loc[all_data.chrom == chrom].copy()
        cpg_bed_sub = cpg_bed.loc[cpg_bed.chrom == chrom].copy()
        
        for _,(chrom, start, end, *_) in cpg_bed_sub.iterrows():
            
            cpgs = all_data_sub.loc[(all_data_sub.start >= start) &
                                    (all_data_sub.start <= end),
                                    all_data_sub.columns.values[3:]].sum(axis=0)
            
            if not cpgs.shape[0]:
                
                continue
            
            count_methylated = cpgs.iloc[range(0, cpgs.shape[0], 2)]
            count_unmethylated = cpgs.iloc[range(1, cpgs.shape[0], 2)]
            
            depth = count_methylated.values + count_unmethylated.values + 1
            methylation_percentage = count_methylated.values / depth
            methylation_percentage[depth < min_depth] = 0
            
            cpg_data = pd.Series(np.concatenate([[chrom, start, end], methylation_percentage]))
            
            methylation_data.append(cpg_data)
            
            all_data_sub = all_data_sub.loc[all_data_sub.start > start,]

    methylation_data = pd.concat(methylation_data, axis=1).T
    
    methylation_data.columns = ['chrom', 'start', 'end'] + list(cov_files.keys())
    
    dtypes = {col : (str if col == 'chrom' else
                     int if col in ['start', 'end'] else
                     float)
              for col in methylation_data.columns}
    
    methylation_data = methylation_data.astype(dtypes)

else:
    
    count_methylated = all_data.iloc[:, range(2, all_data.shape[1], 2)]
    count_unmethylated = all_data.iloc[:, range(3, all_data.shape[1], 2)]
    
    depth = count_methylated.values + count_unmethylated.values + 1
    methylation_percentage = count_methylated.values / depth
    methylation_percentage[depth < min_depth] = 0
    
    methylation_data = pd.DataFrame(np.concatenate([all_data.iloc[:, :2].values,
                                                    methylation_percentage],
                                                   axis=1),
                                    columns = ['chrom', 'start', 'end'] + list(cov_files.keys()))

### Remove elements with 0% methylation in at least 0.25% of the samples

print('Filtering')

thr = 0.25
methylation_data = methylation_data.loc[((methylation_data.iloc[:, 3:] == 0).sum(axis=1) /
                                         (methylation_data.shape[1] - 3)) <= thr,]

### Saving to file

print('Saving to file')

methylation_data.to_csv('methylation_percentage.tsv.gz', sep='\t', header=True, index=False, compression='gzip')
