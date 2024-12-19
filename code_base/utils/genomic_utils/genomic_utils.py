#!/usr/bin/env python3

"""
Collection of useful functions
"""

### ---------------------------------------- ###
### GENOMIC UTILS                            ###
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

def load_gtf(path, desired_biotypes=[], desired_chromosomes=[]):
    
    # Load GTF
    
    gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    gtf_data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Only keep genes

    gtf_data = gtf_data.loc[gtf_data.feature == 'gene', ['seqname', 'start', 'end', 'strand', 'attribute']]
    
    # Get biotype and gene id

    biotypes, gene_ids, gene_symbols = [], [], []
    for _,row in gtf_data.iterrows():
        
        info = row.values[-1]
        
        biotype = re.findall('gene_biotype "\w+";', info)[0]
        biotype = biotype.replace('gene_biotype ', '').replace(';', '').replace('"', '')
        
        biotypes.append(biotype)
        
        gene = re.findall('gene_id "\w+";', info)[0]
        gene = gene.replace('gene_id ', '').replace(';', '').replace('"', '')
        
        gene_ids.append(gene)
        
        if 'gene_name' in info:
            
            gene = info[info.index('gene_name "') + len('gene_name "'):]
            gene = gene[:gene.index('"')]
        
        else:
            
            gene = ''
        
        gene_symbols.append(gene)

    gtf_data['biotype'] = biotypes
    gtf_data['gene_id'] = gene_ids
    gtf_data['gene_symbol'] = gene_symbols
    
    # Filter based on biotype
    
    if len(desired_biotypes):

        gtf_data = gtf_data.loc[gtf_data.biotype.isin(desired_biotypes),]

    # Filter for desired chromosomes

    if len(desired_chromosomes):

        gtf_data = gtf_data.loc[gtf_data.seqname.isin(desired_chromosomes),]
    
    # Remove genes without gene_symbol
    
    #gtf_data = gtf_data.loc[gtf_data['gene_symbol'] != '',]
    
    # Fix dtypes
    
    gtf_data[['start', 'end']] = gtf_data[['start', 'end']].astype(int)
    
    return gtf_data

### ---------------------------------------- ###

def load_vcf(path):
    
    # Get header

    with gzip.open(path,'r') as raw:        
        
        for line in raw:
            
            line = line.decode('utf8')
            
            if line.startswith('#CHROM'):
            
                break
    
    header = line.replace('\n', '').split('\t')
    
    # Load vcf
    
    vcf_data = pd.read_csv(path, sep='\t', header=None, comment='#')

    vcf_data.columns = header
    
    return vcf_data

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

def filter_vcf_with_bed(vcf_data, bed_data):
    
    bed_filter = [bed_data.loc[(bed_data['seqname'] == c) &
                                       (bed_data['start'] <= p) &
                                       (bed_data['end'] >= p),].shape[0] > 0 for _,(c,p) in vcf_data[['#CHROM', 'POS']].iterrows()]
            
    vcf_data = vcf_data.loc[bed_filter,]
    
    return vcf_data
