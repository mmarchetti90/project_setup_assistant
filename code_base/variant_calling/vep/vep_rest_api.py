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
    
    if not cntg.startswith('chr'):
        
        cntg = f'chr{cntg}'
    
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
    
    elif len(r) == len(a):
        
        hgvs_notation = [f'{cntg}:g.{p+n}{ref_base}>{alt_base}' for n,(ref_base,alt_base) in enumerate(zip(r,a)) if ref_base != alt_base and alt_base in 'ACGT']
        
        hgvs_notation = hgvs_notation[0] if len(hgvs_notation) else ''
    
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
        data_parsed['gene_symbol'].append(v['gene_symbol'] if 'gene_symbol' in v.keys() else '')
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
import numpy as np
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

### Summarize variants by most severe effect

# See hierarchy here: https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html

vep_effect_hierarchy = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'feature_elongation',
    'feature_truncation',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_donor_5th_base_variant',
    'splice_region_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'intergenic_variant',
    'sequence_variant'
]

most_severe_effects = []
most_severe_effects_header = ['variant', 'gene_id', 'gene_symbol', 'effect']

for var in np.unique(vars_tsv['variant'].values):
    
    vars_sub = vars_tsv.loc[vars_tsv['variant'] == var,]
    
    affected_genes = np.unique(vars_sub['gene_id'].values)
    
    for gene_id in affected_genes:
        
        vars_sub_sub = vars_sub.loc[vars_sub['gene_id'] == gene_id,]
        
        gene_symbol = vars_sub_sub['gene_symbol'].values[0]
    
        effects = [(c, vep_effect_hierarchy.index(c)) for consequence in vars_sub_sub['consequence_terms'] for c in consequence.split(',')]

        effects.sort(key=lambda c: c[1])
        
        most_severe = effects[0][0]
        
        most_severe_effects.append([var, gene_id, gene_symbol, most_severe])

most_severe_effects = pd.DataFrame(most_severe_effects, columns=most_severe_effects_header)

most_severe_effects.to_csv(f'{output_prefix}_vep_most_severe_effects.tsv', sep='\t', index=False, header=True)
