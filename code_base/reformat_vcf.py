#!/usr/bin/env python3

"""
Class to reformat vcf files in more readable format, selecting desired columns and fields
"""

### ---------------------------------------- ###

class reformat_vcf:
    
    def __init__(self, file_path, desired_columns_and_fields=[]):
        
        # Import file
        data, header, info_fields, columns_names, annotation_fields = self.import_vcf(file_path)
        
        # Subset desired_columns_and_fields
        if desired_columns_and_fields == []:
            
            desired_columns_and_fields = columns_names + info_fields + annotation_fields
        
        else:
            
            desired_columns_and_fields = [dcf for dcf in desired_columns_and_fields if dcf in columns_names + info_fields + annotation_fields]
        
        # Format vcf
        self.format_vcf(data, header, info_fields, columns_names, annotation_fields, desired_columns_and_fields)
        
    ### ------------------------------------ ###
    
    def format_vcf(self, data, header, info_fields, columns_names, annotation_fields, desired_columns_and_fields):
    
        # Init reformatted table
        reformatted_vcf = {}
        
        # Parse data and extract desired info
        for parameter in desired_columns_and_fields:
            
            if parameter in columns_names and parameter != 'INFO':
                
                col_index = columns_names.index(parameter)
                reformatted_vcf[parameter] = [d[col_index] for d in data]
            
            elif parameter in info_fields and parameter != 'ANN':
                
                info_data = []
                
                col_index = columns_names.index('INFO')
                
                for d in data:
                    
                    d_info_data = [field for field in d[col_index].split(';') if field.startswith(f'{parameter}=')]
                    
                    if len(d_info_data) == 1:
                        
                        d_info_data = d_info_data[0].replace(f'{parameter}=', '')
                
                    else:
                        
                        d_info_data = "NA"
                    
                    info_data.append(d_info_data)
                
                reformatted_vcf[parameter] = info_data
            
            elif parameter in annotation_fields:
                
                col_index = columns_names.index('INFO')
                ann_index = annotation_fields.index(parameter)
                
                reformatted_vcf[parameter] = [[field for field in d[col_index].split(';') if field.startswith('ANN=')][0].replace('ANN=', '').split('|')[ann_index] for d in data]
            
            else:
                
                continue
            
        self.reformatted_vcf = pd.DataFrame(reformatted_vcf, columns=desired_columns_and_fields)

    ### ------------------------------------ ###

    def import_vcf(self, file_path):
        
        # Load data
        if file_path.endswith('.gz'):
            
            data = [line.decode() for line in gzip.open(file_path).readlines()]
            
        else:
            
            data = open(file_path).readlines()
        
        # Separate into useful parts
        header = [d for d in data if d.startswith('##')]
        
        info_fields = self.get_info_field_ids([h for h in header if h.startswith('##INFO=')])
        if 'ANN' in info_fields: # snpEff annotations
            
            annotation_fields = [h for h in header if h.startswith('##INFO=')][info_fields.index('ANN')].split('\'')[1].split(' | ')
            
        else:
            
            annotation_fields = []
            
        columns_names = [d for d in data if d.startswith('#') and not d.startswith('##')][0]
        columns_names = columns_names.replace('#', '').replace('\n', '').split('\t')
        
        data = [d.split('\t') for d in data if not d.startswith('#')]
        
        return data, header, info_fields, columns_names, annotation_fields
    
    ### ------------------------------------ ###
    
    def save_reformatted_vcf(self, output_name='reformatted_file.tsv.gz'):
        
        self.reformatted_vcf.to_csv(output_name, sep='\t', header=True, index=False)
    
    ### ------------------------------------ ###
    
    @staticmethod
    def get_info_field_ids(fields):
        
        info_fields = []
        for f in fields:
            
            field_id = f[:f.index(',')].replace('##INFO=<ID=', '')
            info_fields.append(field_id)
            
        return info_fields

### ------------------MAIN------------------ ###

import gzip
import pandas as pd
