#!/usr/bin/env python3

"""
Collection of useful functions
"""

### ---------------------------------------- ###
### GENOMIC DATA INPUT                       ###
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

### ---------------------------------------- ###

def structure_fasta(seq, seq_name='Seq', line_chars=80):
        
    fasta = [f'>{seq_name}']
        
    for i in range(0, len(seq), line_chars):
            
        fasta.append(seq[i : i + line_chars])
        
    fasta = '\n'.join(fasta)
        
    return fasta

### ---------------------------------------- ###

def extract_fasta_kmers(fasta, k=21):
    
    kmers = np.array([])
    
    for contig_name,contig in fasta.items():
        
        new_kmers = np.unique([contig[i : i + k] for i in range(0, len(contig) - k, 1)])
            
        kmers = np.unique(np.concatenate([kmers, new_kmers]))
    
    return kmers

### ---------------------------------------- ###

def load_gtf(path, desired_biotypes=[], desired_chromosomes=[]):
    
    # Load GTF
    
    gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    gtf_data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Only keep genes and exons

    gtf_data = gtf_data.loc[gtf_data.feature.isin(['gene', 'exon']), ['seqname', 'start', 'end', 'strand', 'attribute']]
    
    # Get biotype and gene id

    biotypes, gene_ids, gene_symbols, transcript_ids, exon_ids = [], [], [], [], []
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
        
        if 'transcript_id' in info:
            
            transcript = info[info.index('transcript_id "') + len('transcript_id "'):]
            transcript = transcript[:transcript.index('"')]
        
        else:
            
            transcript = ''
        
        transcript_ids.append(transcript)
        
        if 'exon_id' in info:
            
            exon = info[info.index('exon_id "') + len('exon_id "'):]
            exon = exon[:exon.index('"')]
        
        else:
            
            exon = ''
        
        exon_ids.append(exon)

    gtf_data['biotype'] = biotypes
    gtf_data['gene_id'] = gene_ids
    gtf_data['gene_symbol'] = gene_symbols
    gtf_data['transcript_id'] = transcript_ids
    gtf_data['exon_id'] = exon_ids
    
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
### SEQUENCE ALIGNMENT                       ###
### ---------------------------------------- ###

def gapless_alignment(query, subject, min_overlap=50, matches_reward=1, mismatches_penalty=-1):
    
    # Fix min_overlap if longer than subject
    
    if len(subject) < min_overlap:
        
        min_overlap = len(subject) // 2
    
    # Extend seq1 to allow seq2 to slide on it
    
    ext = len(subject) - min_overlap
    query_extended = ('-' * ext) + query
    
    # Align
    
    seqs_alignments = []
    
    for offset in range(0, len(query_extended) - min_overlap, 1):
        
        matches = sum([matches_reward for a,b in zip(query_extended[offset:], subject) if a == b and a != '-'])
        
        mismatches = sum([mismatches_penalty for a,b in zip(query_extended[offset:], subject) if a != b and (a != '-' or b != '-')])
        
        score = matches + mismatches
        
        seqs_alignments.append((offset, matches, mismatches, score))
    
    # Find best offset
    
    best_offset, max_matches, best_mismatches, best_score = max(seqs_alignments, key=lambda sa: sa[3])
    
    best_offset -= ext
    
    # Define alignment overlap space
    
    """
    N.B. alignment_length = region of overlap between query and seq2 where matches are calculated (e.g. * region)
    N.B. total_length = length of sequence createad by overlapping query and seq2 (e.g. - region)
    
         ATCAGCAGCAGCAGTTGACTGATCGTCAGCTGATCGTACGTGACTGAC
                                        *****************
                                        GATCGTACGTGACTGACACTGACTGTCATCGATCGTACGTCGATGATCGTACGTAC
         ---------------------------------------------------------------------------------------    
    """
    
    alignment_start_query = (best_offset if best_offset >= 0 else 0)
    alignment_start_subject = (0 if best_offset >= 0 else - best_offset)
    
    total_alignment_length = min(len(query) - alignment_start_query, len(subject) - alignment_start_subject)
    
    total_length = max(len(query) + max(0, - best_offset), len(subject) + max(0, best_offset))
    
    query_coverage = total_alignment_length / len(query) # i.e. max overlap, prior to trimming
    
    # Trim alignment at left edge, trying to improve the score

    local_alignment_length = total_alignment_length
    s1, s2 = query[alignment_start_query : alignment_start_query + local_alignment_length], subject[alignment_start_subject : alignment_start_query + local_alignment_length]
    scores = []
    for left_trim in range(local_alignment_length):
        
        mtch = sum([matches_reward for a,b in zip(s1[left_trim:], s2[left_trim:]) if a == b and a != '-'])
        mismtch = sum([mismatches_penalty for a,b in zip(s1[left_trim:], s2[left_trim:]) if a != b and (a != '-' or b != '-')])
        sc = mtch + mismtch
        al_len = local_alignment_length - left_trim
        
        if al_len < min_overlap:
            
            break
        
        scores.append([left_trim, mtch, mismtch, al_len, sc])
        
    scores.sort(key=lambda sc: sc[-1], reverse=True)
    
    best_left_trim = scores[0][0]
    s1, s2 = s1[best_left_trim:], s2[best_left_trim:]
    local_alignment_length = len(s1)

    # Trim alignment at right edge, trying to improve the score

    scores = []
    for right_trim in range(local_alignment_length):
        
        mtch = sum([matches_reward for a,b in zip(s1[:local_alignment_length - right_trim], s2[:local_alignment_length - right_trim]) if a == b and a != '-'])
        mismtch = sum([mismatches_penalty for a,b in zip(s1[:local_alignment_length - right_trim], s2[:local_alignment_length - right_trim]) if a != b and (a != '-' or b != '-')])
        sc = mtch + mismtch
        al_len = local_alignment_length - right_trim
        
        if al_len < min_overlap:
            
            break
        
        scores.append([right_trim, mtch, mismtch, al_len, sc])
        
    scores.sort(key=lambda sc: sc[-1], reverse=True)
    
    best_right_trim = scores[0][0]
    s1, s2 = s1[:len(s1) - best_right_trim], s2[:len(s2) - best_right_trim]
    local_alignment_length = len(s1)
    
    # Finalize score
    
    final_matches = sum([matches_reward for a,b in zip(s1, s2) if a == b and a != '-'])
    final_mismatches = sum([mismatches_penalty for a,b in zip(s1, s2) if a != b and (a != '-' or b != '-')])
    
    final_identity = final_matches / local_alignment_length
    
    final_alignment_score = final_matches * (local_alignment_length / total_length)
    
    return best_offset, final_matches, final_mismatches, final_identity, final_alignment_score, local_alignment_length, total_alignment_length, total_length, query_coverage

### ---------------------------------------- ###

def kmer_alignment(query, subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0):
    
    """
    Local alignment using k-mers matching
    """

    # All uppercase
    
    query, subject = query.upper(), subject.upper()
    
    # Find shared kmers as anchors
    # N.B. using ".index()" when computing match_offsets to make algorithm greedy
    
    query_kmers = [query[i : i + kmer_size] for i in range(0, len(query) - kmer_size + 1, 1)]
    
    subject_kmers = [subject[i : i + kmer_size] for i in range(0, len(subject) - kmer_size + 1, 1)]
    
    match_offsets = np.array([(n, subject_kmers.index(qk), subject_kmers.index(qk) - n)
                              for n,qk in enumerate(query_kmers)
                              if qk in subject_kmers])
    
    # Align
    
    if not len(match_offsets):
        
        return np.array([[], [], []]), (0, 0, 0, 0, 0)
    
    else:
        
        # Align using anchors
        
        alignment = np.repeat('-', 3 * len(query)).reshape((3, len(query)))
        
        alignment[0,] = list(query)
        
        ranges = []
        
        for o in np.unique(match_offsets[:, 2]):
            
            matches_sub = match_offsets[match_offsets[:, 2] == o,].copy()
            
            matches_sub = matches_sub[np.argsort(matches_sub[:, 0])]
            
            if matches_sub.shape[0] < min_anchor_kmers:
                
                continue
            
            q_start, q_end = matches_sub[0, 0], matches_sub[-1, 0]
            
            s_start, s_end = matches_sub[0, 1], matches_sub[-1, 1]
            
            alignment[2, q_start : q_end + kmer_size] = list(subject[s_start : s_end + kmer_size])
        
            ranges.append((q_start, q_end, s_start, s_end))
        
        # Fill in alignment graph
        
        alignment[1, alignment[0,] == alignment[2,]] = '|'
        
        alignment[1, (alignment[0,] != alignment[2,]) & (alignment[2,] != '-')] = ' '
        
        # Alignment stats
        
        alignment_score = (matches_reward * np.sum(alignment[1,] == '|') +
                           mismatches_penalty * np.sum(alignment[1,] == ' ') +
                           gap_penalty * np.sum(alignment[1,] == '-'))
        
        query_coverage = 100 * np.sum(alignment[1,] != '-') / alignment.shape[1]
        
        alignment_length = np.where(alignment[1,] != '-')[0]
        alignment_length = alignment_length[-1] - alignment_length[0] + 1
        
        identity = 100 * np.sum(alignment[1,] == '|') / alignment_length
        
        norm_alignment_score = alignment_score * (alignment_length / alignment.shape[1])
        
        return alignment, (ranges, alignment_length, query_coverage, identity, alignment_score, norm_alignment_score)

### ---------------------------------------- ###

def fill_kmer_alignment_gaps(qry, sbj, graph, stats, kmer_size=5, min_anchor_kmers=5, matches_reward=1, mismatches_penalty=-1, gap_penalty=0):
    
    rng, aln_len, qry_cov, idnt, aln_sc, aln_sc_norm = stats
    
    # Sort ranges
    
    rng.sort(key=lambda r: r[0])
    
    # Fill gaps between anchors
    
    for r1, r2 in zip(rng, rng[1:]):
        
        if r1[1] < r2[0]:
            
            qry_sub = qry[r1[1] : r2[0]]
            
            sbj_sub = sbj[r1[3] : r2[2]]
    
            sub_graph, _ = kmer_alignment(qry_sub, sbj_sub, kmer_size, min_anchor_kmers)
            
            if sub_graph.shape[1] == r2[0] - r1[1]:
            
                graph[:, r1[1] : r2[0]] = sub_graph
    
    # Alignment stats
    
    alignment_score = (matches_reward * np.sum(graph[1,] == '|') +
                       mismatches_penalty * np.sum(graph[1,] == ' ') +
                       gap_penalty * np.sum(graph[1,] == '-'))
    
    query_coverage = 100 * np.sum(graph[1,] != '-') / graph.shape[1]
    
    alignment_length = np.where(graph[1,] != '-')[0]
    alignment_length = alignment_length[-1] - alignment_length[0] + 1
    
    identity = 100 * np.sum(graph[1,] == '|') / alignment_length
    
    norm_alignment_score = alignment_score * (alignment_length / graph.shape[1])
    
    return graph, (rng, alignment_length, query_coverage, identity, alignment_score, norm_alignment_score)

### ---------------------------------------- ###

def structure_kmer_alignment_text(aln, sts, line_chars=60, out_name='alignment.md', save=True):
        
    rng, aln_len, qry_cov, idnt, aln_sc, aln_sc_norm = sts

    # Formatting options for pandoc md to pdf

    header = '\n'.join(['---',
                        'geometry: margin=1.5cm',
                        'papersize: letter',
                        'sansfont:',
                        'fontsize: 12pt',
                        'urlcolor: blue',
                        'toc:',
                        'toc-depth: 4',
                        '---'])
    
    # Format alignment ranges
    
    rng_table = ['## Alignment ranges',
                 '',
                 '| **query_start** | **query_end** | **subject_start** | **subject_end** |',
                 '| :---: | :---: | :---: | :---: |']
    
    rng_table += [f'| {qs} | {qe} | {ss} | {se} |' for qs,qe,ss,se in rng]
    
    rng_table = '\n'.join(rng_table)
    
    # Format stats
    
    stats_table = ['## Stats',
                   '',
                   '| | |',
                   '| :--- | :---: |']
    
    stats_table += [f'| *alignment length* | {aln_len} |',
                    f'| *query coverage* | {qry_cov:.3f}% |',
                    f'| *identity* | {idnt:.3f}% |',
                    f'| *matches* | {(aln[1,] == "|").sum()} |',
                    f'| *mismatches* | {(aln[1,] == " ").sum()} |',
                    f'| *gaps* | {(aln[1,] == "-").sum()} |',
                    f'| *alignment score* | {aln_sc} |',
                    f'| *normalized alignment score* | {aln_sc_norm:.3f} |']
    
    stats_table = '\n'.join(stats_table)
    
    # Format alignment
    
    pos_max_str_len = len(str(aln.shape[1]))
    
    aln_txt = ['## Alignment',
               '',
               f'query 1-{aln.shape[1] + 1}bp',
               '',
               '```latex']
    
    for i in range(0, aln.shape[1], line_chars):
        
        start_pos = str(i + 1)
        start_pos = ' ' * (pos_max_str_len - len(start_pos)) + start_pos + ' '
        
        end_pos = ' ' + str(min(aln.shape[1], i + line_chars + 1))
        
        block = [start_pos + ''.join(aln[0, i : i + line_chars]) + end_pos,
                 ' ' * (pos_max_str_len + 1) + ''.join(aln[1, i : i + line_chars]),
                 ' ' * (pos_max_str_len + 1) + ''.join(aln[2, i : i + line_chars]),
                 '']
        
        aln_txt.extend(block)
    
    aln_txt.append('```')
    
    aln_txt = '\n'.join(aln_txt)
    
    # Finalize outut text, then print
    
    out_txt = '\n\n'.join([header, rng_table, stats_table, aln_txt])
    
    if save:
        
        with open(out_name, 'w') as out_file:
            
            out_file.write(out_txt)
    
    return out_txt

### ---------------------------------------- ###

def reverse_complementary(seq, seq_type='DNA'):
    
    if seq_type == 'DNA':
        
        rev_comp = seq[::-1].lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G')
    
    elif seq_type == 'RNA':
        
        rev_comp = seq[::-1].lower().replace('a', 'U').replace('u', 'A').replace('g', 'C').replace('c', 'G')
    
    else:
        
        rev_comp = ''
        
    return rev_comp

### ---------------------------------------- ###
### VARIANTS                                 ###
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

### ------------------MAIN------------------ ###

try:
    
    import gzip
    import numpy as np
    import pandas as pd
    import re
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()
