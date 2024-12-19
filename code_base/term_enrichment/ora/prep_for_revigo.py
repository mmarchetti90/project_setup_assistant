#!/usr/bin/env python3

"""
This script prepares a GOBP terms table for Revigo from the output of enrichment_analysis_1 and
enrichment_analysis_2
"""

### ---------------------------------------- ###  

def parse_args():

    # GO obo file
    go_file = argv[argv.index('--obo') + 1]

    # Dir containing enrichment files
    enrichment_dir = argv[argv.index('--enrichment_dir') + 1]

    # p threshold
    if '--p_thr' in argv:

        p_thr = float(argv[argv.index('--p_thr') + 1])

    else:

        p_thr = 0.01

    return go_file, enrichment_dir, p_thr

### ---------------------------------------- ###

def parse_go(go_file):
    
    go = {}
    for term in open(go_file).read().split('\n\n'):
        
        if term.startswith('[Term]'):
            
            _, term_id, term_description, namespace, *_ = term.split('\n')
            term_id = term_id.replace('id: ', '')
            term_description = term_description.replace('name: ', '')
            namespace = namespace.replace('namespace: ', '')
            
            if namespace == 'biological_process':
                
                go[term_description] = term_id
    
    return go  

### ------------------MAIN------------------ ###

import pandas as pd

from os import listdir
from sys import argv

### Parse args

go_file, enrichment_dir, p_thr = parse_args()

### Load files

go = parse_go(go_file)

### Create tables

enrichment_analysis_files = [file for file in listdir(enrichment_dir) if file.endswith('.tsv')]

for eaf in enrichment_analysis_files:
    
    file_path = f'{enrichment_dir}/{eaf}'
    
    data = pd.read_csv(file_path, sep='\t')
        
    data = data.loc[data['padj'] < p_thr,]

    data.category = [go[di] for di in data.category.values]

    data = data.loc[:, ['category', 'padj']]
        
    output_name = eaf.replace('.tsv', '_for-revigo.tsv')
    data.to_csv(output_name, sep='\t', header=False, index=False)
