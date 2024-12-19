#!/usr/bin/env python3

"""
Given a vcf file, it runs VEP via its REST API
N.B. Use for human variants only
"""

### ---------------------------------------- ###

def parse_args():
    
    print(argv)
    
    # VCF file
    vcf_path = argv[argv.index('--vcf') + 1]
    
    # Output prefix
    if '--out_prefix' in argv:
        
        output_prefix = argv[argv.index('--out_prefix') + 1]
    
    else:
        
        output_prefix = 'vep_annotation'
    
    # Run SpliceAI?
    if '--spliceai' in argv:
        
        spliceai = True
        
    else:
        
        spliceai = False
    
    return vcf_path, output_prefix, spliceai

### ---------------------------------------- ###

def make_hgvs_notation(r, a, p, cntg):
    
    if len(r) == 1 and len(a) == 1:
        
        # Single base substitution
        hgvs_notation = f'{cntg}:g.{p}{r}>{a}'
    
    elif len(r) < len(a):
        
        # Insertion
        insert_pos = p + len(r) - 1
        insert = a[len(r):]
        insert_len = len(insert)
        hgvs_notation = f'{cntg}:g.{insert_pos}_{insert_pos + 1}ins{insert}'
    
    elif len(r) > len(a):
        
        # Deletion
        deletion_start = p + len(a)
        deletion_end = p + len(r) - 1
        deletion = r[len(a):]
        hgvs_notation = f'{cntg}:g.{deletion_start}_{deletion_end}del{deletion}'
    
    else:
        
        hgvs_notation = ''
    
    return hgvs_notation

### ---------------------------------------- ###

def parse_vars_data_json(dt, spliceai):
    
    data_parsed = {col : []
                   for col in ['variant',
                               'gene_id',
                               'gene_symbol',
                               'transcript_id',
                               'biotype',
                               'impact',
                               'consequence_terms',
                               'DS_AG',
                               'DS_AL',
                               'DS_DG',
                               'DS_DL']}
    
    for v in dt['transcript_consequences']:
        
        data_parsed['variant'].append(dt['input'])
        data_parsed['gene_id'].append(v['gene_id'])
        data_parsed['gene_symbol'].append(v['gene_symbol'])
        data_parsed['transcript_id'].append(v['transcript_id'])
        data_parsed['biotype'].append(v['biotype'])
        data_parsed['impact'].append(v['impact'])
        data_parsed['consequence_terms'].append(','.join(v['consequence_terms']))
        
        if 'spliceai' in v.keys():
            
            data_parsed['DS_AG'].append(v['spliceai']['DS_AG'])
            data_parsed['DS_AL'].append(v['spliceai']['DS_AL'])
            data_parsed['DS_DG'].append(v['spliceai']['DS_DG'])
            data_parsed['DS_DL'].append(v['spliceai']['DS_DL'])
        
        else:
            
            data_parsed['DS_AG'].append('')
            data_parsed['DS_AL'].append('')
            data_parsed['DS_DG'].append('')
            data_parsed['DS_DL'].append('')
    
    data_parsed = pd.DataFrame(data_parsed)
    
    if not spliceai:
        
        data_parsed = data_parsed.loc[:, ~ data_parsed.columns.isin(['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL'])]
    
    return data_parsed

### ------------------MAIN------------------ ###

import json
import pandas as pd
import requests

from sys import argv

### Parse args

vcf_path, output_prefix, spliceai = parse_args()

### Load data

vcf_data = pd.read_csv(vcf_path, sep='\t', header=None, comment='#')

vcf_data = vcf_data.iloc[:, :8]

vcf_data.columns = ['#CHROM',
                    'POS',
                    'ID',
                    'REF',
                    'ALT',
                    'QUAL',
                    'FILTER',
                    'INFO']

### VEP REST calls

vars_json = {}
vars_tsv = []

for _,(contig, pos, _, ref, alt, *_) in vcf_data.iterrows():
    
    for a in alt.split(','):
        
        hgvs_notation = make_hgvs_notation(ref, a, pos, contig)
    
        # Contact VEP REST API
        server = "https://rest.ensembl.org"
        ext = f'/vep/human/hgvs/{hgvs_notation}?'
        
        if spliceai:
            
            ext += 'SpliceAI=1'
        
        rq = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        
        if rq.status_code == 200:
            
            # Get data in Json format, then store in vars_json
            vep_data = rq.json()[0]
            
            vars_json[hgvs_notation] = vep_data
                
            # Parse data for easier read
            vars_data_parsed = parse_vars_data_json(vep_data, spliceai)
            vars_tsv.append(vars_data_parsed)
            
### Dump json with all variants info

with open(f'{output_prefix}_vep.json', 'w') as json_out:

    json.dump(vars_json, json_out, indent=4)

### Save tsv with all variants info

vars_tsv = pd.concat(vars_tsv)

vars_tsv.to_csv(f'{output_prefix}_vep.tsv', sep='\t', index=False, header=True)
