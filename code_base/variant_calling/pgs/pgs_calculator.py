#!/usr/bin/env python3

"""
This script calculates Polygenic Risk Scores from gvcf or vcf files.
N.B. Requires scoring files from https://www.pgscatalog.org/
N.B. Be sure to run gatk HaplotypeCaller with the '-ERC BP_RESOLUTION' option
N.B. Be sure to run gatk GenotypeGVCFs with the '--include-non-variant-sites
true' option
N.B. Alternatively, if the vcf file has only variant sites, a reference fasta
can be provided and scores will be recalculated assuming missing sites to be as
reference
"""

### ---------------------------------------- ###

def parse_args():
    
    # (g)vcf path
    
    vcf_path = argv[argv.index('--vcf') + 1]
    
    # PGS scoring file
    
    pgs_path = argv[argv.index('--pgs') + 1]
    
    # Reference fasta file
    
    if '--fasta' in argv:
    
        fasta_path = argv[argv.index('--fasta') + 1]
    
    else:
        
        fasta_path = 'None'
    
    # Min allele depth
    
    if '--min_depth' in argv:
        
        min_depth = int(argv[argv.index('--min_depth') + 1])
    
    else:
        
        min_depth = 20
    
    return vcf_path, pgs_path, fasta_path, min_depth

### ---------------------------------------- ###

def score_gvcf(path, scoring_dict, min_depth=20):
    
    header = []
    
    if path.endswith('.gz'):
        
        raw_input = gzip.open(path, 'r')
    
    else:
        
        raw_input = open(path, 'r')
    
    for line in raw_input:
        
        if type(line) == bytes:
            
            line = line.decode()
        
        line = line.strip()
    
        line = line.replace('\n', '')
        
        if line.startswith('##'):
                
            header.append(line)
        
        elif line.startswith('#CHROM'):
            
            columns = line.split('\t')
            
            pgs_scores = {c : 0 for c in columns[9:]}
            
            found_alleles = {c : set() for c in columns[9:]}
            
        elif len(line):
            
            var_dict = {c : v for c,v in zip(columns, line.split('\t'))}
            
            # Only keep reference and SNPs
            
            alleles = get_alleles(var_dict['REF'], var_dict['ALT'].split(','), var_dict['POS'], var_dict['#CHROM'])
            
            # Skip if no allele is in the PGS scoring dictionary
            
            alleles_in_scoring_dict = [a for allele in alleles.values() for a in allele if a in scoring_dict.keys()]
            
            if len(alleles_in_scoring_dict) > 0:
                
                gt_col = var_dict['FORMAT'].split(':').index('GT')
                
                if 'AD' in var_dict['FORMAT'].split(':'):
                
                    ad_col = var_dict['FORMAT'].split(':').index('AD')
                
                else:
                    
                    ad_col = -1
                
                for sample in pgs_scores.keys():
                    
                    sample_gt = var_dict[sample].split(':')[gt_col].replace('|', '/')
                    
                    if ad_col >= 0:
                    
                        sample_ad = var_dict[sample].split(':')[ad_col].split(',')
                    
                    else:
                        
                        sample_ad = [min_depth + 1, min_depth + 1]
                    
                    for sgt in sample_gt.split('/'):
                        
                        # Check if allele is valid
                        
                        if sgt not in alleles.keys():
                            
                            continue
                        
                        # Check if allele has enough coverage
                        
                        try:
                            
                            allele_ad = int(sample_ad[int(sgt)])
                            
                        except:
                            
                            allele_ad = 0
                        
                        if allele_ad < min_depth:
                            
                            continue
                        
                        # Update sample score
                        
                        for allele in alleles[sgt]:
                            
                            if allele in scoring_dict.keys():
                                
                                pgs_scores[sample] += scoring_dict[allele]
                                
                                found_alleles[sample].add(allele)
        
        else:
            
            pass
    
    return pgs_scores, found_alleles

### ---------------------------------------- ###

def get_alleles(ref, alts, pos, cntg):
    
    pos = int(pos)

    alleles = {}
    
    # Reference
    
    allele_id = '0'
    
    alleles[allele_id] = [encode_allele(cntg, pos + i, r) for i,r in enumerate(ref) if r in 'ACGT']
    
    # Alternatives
    
    for i,alt in enumerate(alts):
        
        allele_id = str(i + 1)
        
        if alt == '<*>':
            
            continue
        
        if len(alt) == 1:
            
            alleles[allele_id] = [encode_allele(cntg, pos, alt)]
        
        elif len(ref) < len(alt):
            
            # Insertion, skipped
            
            alleles[allele_id] = ['']
        
        elif len(ref) > len(alt):
            
            # Deletion, skipped
            
            alleles[allele_id] = ['']
        
        elif len(ref) == len(alt):
            
            alleles[allele_id] = [encode_allele(cntg, pos + i, a) for i,a in enumerate(alt) if a in 'ACGT']
        
        else:
            
            alleles[allele_id] = ['']
    
    return alleles

### ---------------------------------------- ###

def encode_allele(cntg, pos, allele):
    
    if not cntg.startswith('chr'):
        
        cntg = f'chr{cntg}'
    
    encoded = f'{cntg}:{pos}:{allele}'
    
    return encoded

### ---------------------------------------- ###

def decode_allele(encoded):
    
    cntg, pos, allele = encoded.split(':')
    
    pos = int(pos)
    
    return cntg, pos, allele

### ---------------------------------------- ###

def load_fasta(path):
    
    if path.endswith('.gz'):
        
        chromosomes = gzip.open(path, 'r').read().decode().split('>')
    
    else:
        
        chromosomes = open(path, 'r').read().split('>')
    
    fasta = {}
    
    for chrom in chromosomes:
        
        if not len(chrom):
            
            continue
        
        chrom = chrom.split('\n')
        
        chrom_name = chrom[0].split(' ')[0]
        
        chrom_seq = ''.join(chrom[1:])
        
        fasta[chrom_name] = chrom_seq
    
    return fasta

### ------------------MAIN------------------ ###

import gzip
import pandas as pd

from sys import argv

### Parse args

vcf_path, pgs_path, fasta_path, min_depth = parse_args()

### Load PGS scoring file, then create dict with hgvs-like variants and effects as keys and values, respectively

pgs = pd.read_csv(pgs_path, sep='\t', comment='#', dtype=str)

pgs_dict = {encode_allele(row['chr_name'], row['chr_position'], row['effect_allele']) : float(row['effect_weight']) for _,row in pgs.iterrows()}

### Parse (g)vcf and calculate scores

pgs_scores, found_alleles = score_gvcf(vcf_path, pgs_dict, min_depth)

### If a reference was provided, add information about non-variant sites

if fasta_path != 'None':
    
    pgs_scores_backup = pgs_scores.copy()
    
    # Load contigs
    
    contigs = load_fasta(fasta_path)
    
    # Compute list of alleles that match the reference fasta
    
    reference_pgs_alleles = []
    
    for allele in pgs_dict.keys():
        
        cntg, pos, pgs_allele = decode_allele(allele)
        
        if cntg in contigs.keys():
            
            if pos > len(contigs[cntg]): # Sanity check
                
                continue
            
            ref_allele = contigs[cntg][pos - 1] # VCF are 1 based but contigs is 0 based
            
        elif cntg.replace('chr', '') in contigs.keys():
            
            cntg = cntg.replace('chr', '')
            
            if pos > len(contigs[cntg]): # Sanity check
                
                continue
            
            ref_allele = contigs[cntg][pos - 1] # VCF are 1 based but contigs is 0 based
        
        else:
            
            continue
        
        if ref_allele == pgs_allele:
            
            reference_pgs_alleles.append(allele)
    
    # Update scores
    
    for sample in pgs_scores.keys():
        
        sample_score, sample_alleles = pgs_scores[sample], found_alleles[sample]
        
        missing_alleles = [allele for allele in reference_pgs_alleles if allele not in sample_alleles]
        
        sample_score += sum([2 * pgs_dict[ma] for ma in missing_alleles])
        
        pgs_scores[sample] = sample_score

### Reformat as data frames

pgs_scores = pd.DataFrame.from_dict(pgs_scores, orient='index').reset_index(drop=False)
pgs_scores.columns = ['sample', 'pgs']

found_alleles = pd.DataFrame.from_dict(found_alleles, orient='index').T

### Save data to file

pgs_scores.to_csv('pgs_scores.tsv', sep='\t', index=False)

found_alleles.to_csv('found_alleles.txt', sep='\t', index=False)
