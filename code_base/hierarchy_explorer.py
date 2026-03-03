#!/usr/bin/env python3

### ---------------------------------------- ###

class explore_hierarchy:
    
    """
    Class to quickly parse a hierarchy database in the form of a tab-delimited file with two
    columns, "parent" and "child", in that order
    (N.B. first row must be "parent\tchild")
    """
    
    def __init__(self, hierarchy_file):
        
        self.hierarchy = pd.read_csv(hierarchy_file, sep='\t')

    ### ------------------------------------ ###
    ### UTILITIES                            ###
    ### ------------------------------------ ###
    
    def get_children_terms(self, c):
        
        children = set(self.hierarchy.loc[self.hierarchy.parent.isin(c), 'child'].to_list())
        
        children.update(c)
    
        if children == c:
            
            return c
        
        else:
            
            return self.get_children_terms(children)
    
    ### ------------------------------------ ###
    
    def get_parent_terms(self, p):
        
        parents = set(self.hierarchy.loc[self.hierarchy.child.isin(p), 'parent'].to_list())
        
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
            
            children = set(self.hierarchy.loc[self.hierarchy.parent.isin(children), 'child'].to_list())
        
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
