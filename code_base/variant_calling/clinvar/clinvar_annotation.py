#!/usr/bin/env python3

"""
Interpolating a target vcf to clinvar's vcf

N.B. Download Clinvar's vcf from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
     or set '--clinvar_vcf' to 'download'
"""

### ---------------------------------------- ###

def parse_args():
    
    print(argv)
    
    # Target VCF file
    
    target_vcf_path = argv[argv.index('--target_vcf') + 1]
    
    target_vcf_data = pd.read_csv(target_vcf_path, sep='\t', header=None, comment='#', dtype=str)
    
    target_vcf_data = target_vcf_data.iloc[:, :8]

    target_vcf_data.columns = [
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    ]
    
    # Clinvar VCF file
    
    clinvar_vcf_path = argv[argv.index('--clinvar_vcf') + 1]
    
    if clinvar_vcf_path == 'download':
        
        url = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
        
        clinvar_vcf_data = pd.read_csv(url, sep='\t', header=None, comment='#', dtype=str)
    
    else:
        
        clinvar_vcf_data = pd.read_csv(clinvar_vcf_path, sep='\t', header=None, comment='#', dtype=str)
    
    clinvar_vcf_data = clinvar_vcf_data.iloc[:, :8]

    clinvar_vcf_data.columns = [
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    ]
    
    return target_vcf_data, clinvar_vcf_data

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

### ------------------MAIN------------------ ###

import pandas as pd

from sys import argv

### Load data

target_vcf_data, clinvar_vcf_data = parse_args()

### Generate hgvs names for Clinvar's vcf

clinvar_vcf_data.loc[:, 'hgvs'] = [make_hgvs_notation(ref, alt, int(pos), contig) for _,(contig, pos, _, ref, alt, *_) in clinvar_vcf_data.iterrows()]

### Parse target variants

vars_annotation = []
vars_annotation_header = ['hgvs', 'contig', 'pos', 'ref', 'alt', 'affected_genes', 'clinical_significance', 'associated_disease']

for _,(contig, pos, _, ref, alt, *_) in target_vcf_data.iterrows():
    
    for a in alt.split(','):
        
        hgvs_notation = make_hgvs_notation(ref, a, int(pos), contig)
        
        if hgvs_notation in clinvar_vcf_data['hgvs'].values:
            
            clinvar_sub = clinvar_vcf_data.loc[clinvar_vcf_data['hgvs'] == hgvs_notation, 'INFO'].values
            
            affected_genes = ';'.join([i.replace('GENEINFO=', '') for info in clinvar_sub for i in info.split(';') if i.startswith('GENEINFO=')])
            
            clinical_significance = ';'.join([i.replace('CLNSIG=', '') for info in clinvar_sub for i in info.split(';') if i.startswith('CLNSIG=')])
            
            associated_disease = ';'.join([i.replace('CLNDN=', '') for info in clinvar_sub for i in info.split(';') if i.startswith('CLNDN=')])
        
        else:
            
            affected_genes = 'unknown'
            
            clinical_significance = 'not_in_clinvar'
            
            associated_disease = 'unknown'

        vars_annotation.append([hgvs_notation, contig, pos, ref, alt, affected_genes, clinical_significance, associated_disease])

vars_annotation = pd.DataFrame(vars_annotation, columns=vars_annotation_header)

vars_annotation.to_csv('clinvar_annotation.tsv', sep='\t', index=False, header=True)
