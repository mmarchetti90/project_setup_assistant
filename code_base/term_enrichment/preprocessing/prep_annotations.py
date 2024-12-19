#!/usr/bin/env python3

"""
This script prepares a term to gene table for GOBP terms
"""

### ---------------------------------------- ###

def parse_args():
    
    # GO file
    obo_file = argv[argv.index('--obo') + 1]
    
    # Annotation file
    gaf_file = argv[argv.index('--gaf') + 1]
    
    # GO file
    namespace = argv[argv.index('--namespace') + 1]
    
    return obo_file, gaf_file, namespace

### ---------------------------------------- ###

def parse_go(go_file, target_namespace):

    # This creates a tree structure for deep search
    
    # Read obo file and create hierarchy
    id_to_description = {}
    go_hierarchy = {}
    for term in open(go_file).read().split('\n\n'):
        
        if term.startswith('[Term]'):
            
            # Parse term info
            _, term_id, term_description, namespace, *info = term.split('\n')
            term_id = term_id.replace('id: ', '')
            term_description = term_description.replace('name: ', '')
            namespace = namespace.replace('namespace: ', '')
            
            # Skip if term is obsolete
            if 'is_obsolete: true' in info:
                
                continue
            
            # Skip if wrong namespace
            if namespace != target_namespace:
                
                continue
            
            # Store description
            id_to_description[term_id] = term_description
            
            # Store root hierarchy
            if term_id not in go_hierarchy.keys():
                
                go_hierarchy[term_id] = []
            
            # Find parent terms, if they exist, then add the term as a child
            parents = ([i.split(' ')[1] for i in info if i.startswith('is_a:')] +
                       [i.split(' ')[2] for i in info if i.startswith('relationship:')])
            
            if len(parents):
                
                for p in parents:
                    
                    try:
                        
                        go_hierarchy[p].append(term_id)
                    
                    except:
                        
                        go_hierarchy[p] = [term_id]
    
    # Remove duplicates
    go_hierarchy = {term : list(set(hierarchy)) for term,hierarchy in go_hierarchy.items()}
    
    # Convert hierarchy to a faster file format
    go_hierarchy_fast = []
    for term in go_hierarchy.keys():
        
        parents = set([p for p,c in go_hierarchy.items() if term in c])
        
        for p in parents:
            
            go_hierarchy_fast.append((p, term))

    go_hierarchy_fast = pd.DataFrame(go_hierarchy_fast, columns=['parent', 'child'])
    
    return id_to_description, go_hierarchy_fast

### ---------------------------------------- ###

def create_term2gene(id_to_description, go_hierarchy, annotation):
    
    # Subset annotation (keep only entries with terms in go_hierarchy, then remove unwanted columns)
    annotation = annotation.loc[annotation.iloc[:, 4].isin(go_hierarchy.values.flatten()),]
    annotation = pd.DataFrame({'gs_name' : annotation.iloc[:, 4].to_list(),
                               'gene_name' : annotation.iloc[:, 2].to_list()})
    
    # Init term2gene
    term2gene = {'go_id' : [],
                 'gs_name' : [],
                 'gene_name' : []}
    
    # Populate term2gene
    for term in set(go_hierarchy.values.flatten()):
        
        # Get term description
        description = id_to_description[term]
        
        # Extract all children terms
        all_terms = get_children_terms({term}, go_hierarchy)
        
        # Get genes associated with term or its children
        genes = list(set(annotation.loc[annotation.gs_name.isin(all_terms), 'gene_name']))
        
        # Add data to term2gene
        term2gene['go_id'].extend([term for _ in range(len(genes))])
        term2gene['gs_name'].extend([description for _ in range(len(genes))])
        term2gene['gene_name'].extend(genes)
        
    term2gene = pd.DataFrame(term2gene)
    
    term2gene.sort_values(by='gs_name', inplace=True)
    
    term2gene.drop_duplicates(inplace=True)
    
    return term2gene

### ---------------------------------------- ###

def get_children_terms(c, parent_child):
    
    children = set(parent_child.loc[parent_child.parent.isin(c), 'child'].to_list())
    
    children.update(c)

    if children == c:
        
        return c
    
    else:
        
        return get_children_terms(children, parent_child)

### ------------------MAIN------------------ ###

import pandas as pd

from sys import argv

### Read args

go_file, annotation_file, namespace = parse_args()

### Load files

id_to_description, go_hierarchy = parse_go(go_file, namespace)
annotation = pd.read_csv(annotation_file, sep='\t', comment='!', header=None, dtype=str)

### Create term2gene file

term2gene = create_term2gene(id_to_description, go_hierarchy, annotation)
term2gene = term2gene.loc[:, term2gene.columns != 'go_id']
term2gene.to_csv("term2gene.tsv", sep='\t', index=False, header=True)
