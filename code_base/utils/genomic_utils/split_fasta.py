#!/usr/bin/env python3

"""
Converts a multi-contig fasta to individual fasta files
"""

### ---------------------------------------- ###

def load_fasta(path):
    
    fasta = {}
    
    for chrom in open(path).read().split('>'):
        
        if not len(chrom):
            
            continue
        
        chrom = chrom.split('\n')
        
        chrom_name = chrom[0].split(' ')[0]
        
        chrom_seq = ''.join(chrom[1:])
        
        fasta[chrom_name] = chrom_seq
    
    return fasta

### ---------------------------------------- ###

def structure_fasta(seq, seq_name='Seq', line_chars=80):
        
    fasta = [f'>{seq_name}']
        
    for i in range(0, len(seq), line_chars):
            
        fasta.append(seq[i : i + line_chars])
        
    fasta = '\n'.join(fasta)
        
    return fasta

### ------------------MAIN------------------ ###

from sys import argv

fasta_file = argv[argv.index('--fasta') + 1]

fasta = load_fasta(fasta_file)

print('contig\tsize(bp)')

for contig_name, contig_seq in fasta.items():

    print(f'{contig_name}\t{len(contig_seq)}')
    
    fasta_out_txt = structure_fasta(contig_seq, contig_name, 80)

    with open(f'{contig_name}.fasta', 'w') as fasta_out_file:
        
        fasta_out_file.write(fasta_out_txt)
