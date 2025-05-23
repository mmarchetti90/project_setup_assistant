#!/usr/bin/env python3"""This script checks for strandness of gene counts from STAR --quantMode"""### ---------------------------------------- ###def analyse_strandness(f):        data = pd.read_csv(f, sep='\t', header=None, index_col=0)    data.columns = ['unstranded', 'firststrand', 'secondstrand']        _, f1_nofeature, f2_nofeature = data.loc['N_noFeature',]    _, f1_withfeature, f2_withfeature = data.iloc[4:].sum(axis=0)        f1_percent, f2_percent = 100 * f1_withfeature / (f1_withfeature + f1_nofeature), 100 * f2_withfeature / (f2_withfeature + f2_nofeature)        f1_p = np.exp(f1_percent) / (np.exp(f1_percent) + np.exp(f2_percent))    f2_p = 1 - f1_p        if f1_p > f2_p and f1_p > 0.9:                strandness = 'firststrand'            elif f2_p > f1_p and f2_p > 0.9:                strandness = 'secondstrand'            else:                'unclear'        return f1_percent, f2_percent, f1_p, f2_p, strandness### ------------------MAIN------------------ ###import numpy as npimport pandas as pdfrom os import listdirfrom sys import argvfile_suffix = argv[argv.index('--file_suffix') + 1]### Get file listfiles = [file for file in listdir() if file.endswith(file_suffix)]files.sort()### Parse filesstrand_info = {'Sample' : [],               'FirststrandWithFeature%' : [],               'SecondstrandWithFeature%' : [],               'FirststrandP' : [],               'SecondstrandP' : [],               'Strandness' : []}for file in files:        sample = file.replace(file_suffix, '')        f1_percent, f2_percent, f1_p, f2_p, strandness = analyse_strandness(file)        strand_info['Sample'].append(sample)    strand_info['FirststrandWithFeature%'].append(f1_percent)    strand_info['SecondstrandWithFeature%'].append(f2_percent)    strand_info['FirststrandP'].append(f1_p)    strand_info['SecondstrandP'].append(f2_p)    strand_info['Strandness'].append(strandness)    strand_info = pd.DataFrame(strand_info)strand_info.to_csv('strand_info.tsv', sep='\t', index=False, header=True)