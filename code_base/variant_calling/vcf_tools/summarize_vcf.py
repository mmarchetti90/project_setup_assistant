#!/usr/bin/env python3

"""
This script parses the variant calling results, removes variants from control samples and outputs a report
"""

### ---------------------------------------- ###

def parse_args():
    
    # VCF parameters to output
    if '--desired_parameters' in argv:
        
        desired_columns_and_fields = argv[argv.index('--desired_parameters') + 1].split(';')
        
        for necessary_field in ['CHROM', 'POS', 'REF', 'ALT']:
            
            if necessary_field not in desired_columns_and_fields:
                
                desired_columns_and_fields = [necessary_field] + desired_columns_and_fields
        
    else:
        
        desired_columns_and_fields = []
    
    # Sample file
    sample_file_path = argv[argv.index('--sample') + 1]
    sample_name = basename(sample_file_path).replace('.vcf', '').replace('.gz', '')
    sample_vcf = rv.reformat_vcf(sample_file_path, desired_columns_and_fields).reformatted_vcf
    
    # Control file(s)
    if '--control' in argv:
        
        control_file_paths = argv[argv.index('--control') + 1].split(',')
        
    else:
        
        control_file_paths = []
    
    control_vcf = [rv.reformat_vcf(path, desired_columns_and_fields).reformatted_vcf for path in control_file_paths]
    
    return sample_name, sample_vcf, control_vcf

### ---------------------------------------- ###

def clean_vcf(sample_vcf, control_vcf):
    
    # Removes variants found in control samples
    
    required_fields = ['CHROM', 'POS', 'REF', 'ALT']
    
    # Make list of variants in sample
    sample_variants = sample_vcf.loc[:, required_fields].apply(lambda x: "_".join(x), axis =1).to_list()
    
    # Make list of variants in controls
    ctrl_variants = []
    
    for ctrl in control_vcf:
        
        ctrl_variants.extend(ctrl.loc[:, required_fields].apply(lambda x: "_".join(x), axis =1).to_list())
    
    # Filter sample variants
    variants_filter = [v not in ctrl_variants for v in sample_variants]
    sample_vcf = sample_vcf.loc[variants_filter,]
    
    return sample_vcf

### ------------------MAIN------------------ ###

import reformat_vcf as rv

from os.path import basename
from sys import argv

# Parse arguments
sample_name, sample_vcf, control_vcf = parse_args()

# Remove variants found in control samples
sample_vcf = clean_vcf(sample_vcf, control_vcf)

# Export data
if len(control_vcf):

    sample_vcf.to_csv(f'{sample_name}_cleaned.tsv.gz', sep='\t', header=True, index=False)

else:

    sample_vcf.to_csv(f'{sample_name}.tsv.gz', sep='\t', header=True, index=False)
