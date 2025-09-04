#!/usr/bin/env python3

### ---------------------------------------- ###

class explore_taxonomy:
    
    """
    Class to parse a taxonomy hierarchy database prepared using prepare_taxonomy_db.py
    
    Parameters:
    taxid_hierarchy_file
        Path to file containing hierarchy information stored as tab-delimited text with two
        columns, "parent" and "child", in that order
        (N.B. first row must be "parent\tchild")
    taxid_to_name_file
        Path to file containing a taxid to name conversion stored as tab-delimited text with two
        columns, "taxid" and "name"
        (N.B. column names are omitted)
    taxid_to_rank_file
        Path to file containing a taxid to rank conversion stored as tab-delimited text with two
        columns, "taxid" and "rank"
        (N.B. column names are omitted)
    """
    
    def __init__(self, taxid_hierarchy_file, taxid_to_name_file='', taxid_to_rank_file=''):
        
        self.taxid_hierarchy = pd.read_csv(taxid_hierarchy_file, sep='\t')
        
        try:
        
            self.taxid_to_name = {int(row.split('\t')[0]) :row.split('\t')[1]
                                  for row in open(taxid_to_name_file).read().split('\n')
                                  if len(row)}
        
        except:
            
            self.taxid_to_name = {}
        
        try:
            
            self.taxid_to_rank = {int(row.split('\t')[0]) :row.split('\t')[1]
                                  for row in open(taxid_to_rank_file).read().split('\n')
                                  if len(row)}

        except:
            
            self.taxid_to_rank = {}

    ### ------------------------------------ ###
    ### UTILITIES                            ###
    ### ------------------------------------ ###
    
    def get_children_terms(self, c):
        
        children = set(self.taxid_hierarchy.loc[self.taxid_hierarchy.parent.isin(c), 'child'].to_list())
        
        children.update(c)
    
        if children == c:
            
            return c
        
        else:
            
            return self.get_children_terms(children)
    
    ### ------------------------------------ ###
    
    def get_parent_terms(self, p):
        
        parents = set(self.taxid_hierarchy.loc[self.taxid_hierarchy.child.isin(p), 'parent'].to_list())
        
        parents.update(p)
    
        if parents == p:
            
            return p
        
        else:
            
            return self.get_parent_terms(parents)
    
    ### ------------------------------------ ###
    
    def get_parent_child_distance(self, p, c):
        
        children = [p]
        
        n = 0
        
        while c not in children:
            
            n += 1
            
            children = set(self.taxid_hierarchy.loc[self.taxid_hierarchy.parent.isin(children), 'child'].to_list())
        
        return n
    
    ### ------------------------------------ ###
    
    def find_closest_common_ancestor(self, element_1, element_2):
        
        # Get ancestors in common
        
        parents_1 = self.get_parent_terms([element_1])
        
        parents_2 = self.get_parent_terms([element_2])
        
        common_ancestors = [p1 for p1 in parents_1 if p1 in parents_2 and p1 != '']
        
        # Find shortest distance
        
        elements_distance = []
        
        for ca in common_ancestors:
            
            distance_1 = self.get_parent_child_distance(ca, element_1)
            
            distance_2 = self.get_parent_child_distance(ca, element_2)
            
            elements_distance.append((ca, distance_1 + distance_2 + 1))
        
        elements_distance.sort(key=lambda ed: ed[1])
        
        # Closest ancestor
        
        closest_ancestor = elements_distance[0]
        
        return closest_ancestor

### ---------------------------------------- ###

import pandas as pd
