#!/usr/bin/env python3

"""
This script reads in a gene list and extract DNA sequences around their promoter
"""

### ---------------------------------------- ###

def parseArgs():
    
    # import gene list
    gene_list_path = argv[argv.index("--gene_list") + 1]
    gene_list = open(gene_list_path).read().split('\n')
    
    # Genome fasta
    genome_fasta_path = argv[argv.index("--genome_fasta") + 1]
    if genome_fasta_path.endswith('.gz'):
        
        genome_fasta = gzip.open(genome_fasta_path).read().decode('utf8').split('>')
    
    else:
        
        genome_fasta = open(genome_fasta_path).read().split('>')
    
    genome_fasta = parseFasta(genome_fasta)
    
    # Genome annotation
    genome_annotation_path = argv[argv.index("--genome_annotation") + 1]
    if genome_annotation_path.endswith('.gz'):
        
        genome_annotation = gzip.open(genome_annotation_path).read().decode('utf8').split('\n')
    
    else:
        
        genome_annotation = open(genome_annotation_path).read().split('\n')
    
    genome_annotation = parseAnnotation(genome_annotation)
    
    # Upstream window
    if '--upstream' in argv:
        
        upstream = int(argv[argv.index("--upstream") + 1])
        
    else:
        
        upstream = 2000
    
    # Downstream window
    if '--downstream' in argv:
        
        downstream = int(argv[argv.index("--downstream") + 1])
        
    else:
        
        downstream = 500
    
    return gene_list, genome_fasta, genome_annotation, upstream, downstream

### ---------------------------------------- ###

def parseFasta(genome_fasta):
    
    parsed = {}
    for gf in genome_fasta:
        
        if not len(gf):
            
            continue
        
        chromosome = gf.split('\n')[0].split(' ')[0]
        sequence = ''.join(gf.split('\n')[1:])
        
        parsed[chromosome] = sequence
    
    return parsed

### ---------------------------------------- ###

def parseAnnotation(genome_annotation):
    
    parsed = {}
    for ga in genome_annotation:
        
        if not len(ga):
            
            continue
        
        chromosome, _, element_type, start, stop, _, strand, _, info = ga.split('\t')
        
        if element_type != 'gene':
            
            continue
        
        gene_id = extractInfo(info, 'gene_id "')
        gene_symbol = extractInfo(info, 'gene_symbol "')
        
        parsed[gene_id] = {'Chr' : chromosome,
                           'Start' : start,
                           'Stop' : stop,
                           'Strand' : strand}
        
        parsed[gene_symbol] = {'Chr' : chromosome,
                               'Start' : start,
                               'Stop' : stop,
                               'Strand' : strand}
    
    return parsed

### ---------------------------------------- ###

def extractInfo(text, field):
    
    requested_info = text[text.index(field) + len(field) : text[text.index(field) : ].index('";') + text.index(field)]
    
    return requested_info

### ---------------------------------------- ###

def extractSequences(gene_list, genome_fasta, genome_annotation, upstream, downstream):
    
    sequences = []
    for (i, gene) in enumerate(gene_list):
        
        if gene not in genome_annotation.keys():
            
            continue
        
        chromosome, start, stop, strand = genome_annotation[gene].values()
        
        if chromosome not in genome_fasta.keys():
            
            continue
        
        if strand == '+':
            
            sequence = genome_fasta[chromosome][max(0, int(start) - int(upstream)) : int(start) + int(downstream)]
        
        else:
            
            sequence = genome_fasta[chromosome][max(0, int(start) - int(downstream)) : int(start) + int(upstream)]
        
        sequences.append(f'>{gene}')
        sequences.append(sequence)
        sequences.append('\n')
        
    sequences = '\n'.join(sequences)
    
    with open('promoter_sequences.fa', 'w') as fasta_out:
        
        fasta_out.write(sequences)

### ------------------MAIN------------------ ###

import gzip

from sys import argv

# Import data
gene_list, genome_fasta, genome_annotation, upstream, downstream = parseArgs()

# Extract sequences
extractSequences(gene_list, genome_fasta, genome_annotation, upstream, downstream)
