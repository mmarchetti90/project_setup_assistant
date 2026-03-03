#!/usr/bin/env python3

### ---------------------------------------- ###

class cluster_terms:
    
    """
    Class to clusters hierarchical terms (e.g. HPO, Gene Ontology, Reactome, etc) based on their distance
    Only terms in the same branch are clustered
    Clusters are based on graph leaves (i.e. it used terms with no childrens among the other terms)

    N.B. Algorithm is stringent

    Arguments
    ---------
    terms_list : list[str]
        List of terms strings
    parent_child : pandas frame
        Pandas data frame of term hierarchy
        Must have two columns: parent and child
    max_distance
        Maximum distance between terms to allow clustering
        (Default = 3)
    cross_branches : bool
        Set to True to allow distance calculation to cross branches through a common ancestor
        (Default = True)
    output_prefix : str
        Prefix for output files
        (Default = 'terms')
    """

    def __init__(self, terms_list, parent_child, max_distance=3, cross_branches=True, output_prefix='terms'):
        
        self.terms_list = terms_list
        
        self.parent_child = parent_child
        
        self.max_distance = max_distance
        
        self.cross_branches = cross_branches
        
        self.output_prefix = output_prefix

    ### ------------------------------------ ###
    ### CLUSTERING                           ###
    ### ------------------------------------ ###
    
    def cluster(self):
        
        # Compute pairwise distances
        
        print('Computing pairwise distances.')

        self.term_distance()

        # Cluster
        
        print('Clustering terms.')

        self.cluster_terms()

        # Simplifying terms
        
        print('Simplifying terms')

        self.simplified_terms = list(self.clusters_parents.values())
        
        # Saving data

        self.save_data()

    ### ------------------------------------ ###
    
    def term_distance(self):
        
        terms = self.terms_list
        
        impossible_distance = -1
        
        distances = np.zeros((len(terms), len(terms)))
        
        for t1 in terms:
            
            idx1 = terms.index(t1)
            
            t1_parents = self.get_parent_terms([t1])
            
            t1_children = self.get_children_terms([t1])
            
            for t2 in terms[terms.index(t1) + 1:]:
                
                idx2 = terms.index(t2)
                
                if t2 in t1_parents:
                    
                    dist = self.get_parent_child_distance(t2, t1)
                
                elif t2 in t1_children:
                    
                    dist = self.get_parent_child_distance(t1, t2)
                
                else:
                    
                    if self.cross_branches:
                    
                        _, dist = self.find_closest_common_ancestor(t1, t2)
                    
                    else:
                        
                        dist = impossible_distance
                
                distances[idx1, idx2] = dist
                
                distances[idx2, idx1] = dist
                
        self.pairwise_distances = pd.DataFrame(distances, index=terms, columns=terms)

    ### ------------------------------------ ###
    
    def cluster_terms(self):
        
        terms = self.terms_list
        
        log = []
        
        # Find childrenless terms (i.e. leaves)
        # N.B. a node with childrens is still considered a leaf if none of the children is in terms
        
        leaves = []
        
        for h in terms:
            
            childrens = [c for c in self.get_children_terms([h]) if c in terms]
            
            if len(childrens) == 1:
                
                leaves.append(h)
        
        original_leaves_num = len(leaves)
        
        log.append(f'Found {original_leaves_num} initial leaves')
        
        # Temporarily remove terms shared by leaves
        
        sharing = {}
        
        for l in leaves:
            
            parents = self.get_parent_terms([l])
            
            for p in parents:
                
                if p not in terms:
                    
                    continue
                
                if p in sharing.keys():
                    
                    sharing[p] += 1
                
                else:
                    
                    sharing[p] = 1
        
        to_be_removed = [s for s,sn in sharing.items() if sn > 1]
        
        terms = [h for h in terms if h not in to_be_removed]
        
        log.append(f'Temporarily removed {len(to_be_removed)} terms shared between leaves ancestors')
        
        # Init clusters
        
        clusters = {h : (leaves.index(h) if h in leaves
                         else -1)
                    for h in terms}
        
        clusters_parents = {cl : h for h,cl in clusters.items() if cl != -1}
        
        # Cluster
        
        toggle = True
        
        while toggle:
            
            old_clusters = clusters.copy()
            
            for l in leaves:
                
                # Extract cluster elements
                
                l_clust = clusters[l]
                
                cluster_elements = [h for h,cl in clusters.items() if cl == l_clust]
                
                cluster_parent = clusters_parents[l_clust]
                
                # Find closest non-cluster element which is a parent of cluster_parent
                
                closest_term = [(p, self.pairwise_distances.loc[p, cluster_parent]) for p in self.get_parent_terms([cluster_parent]) if p in terms and p not in cluster_elements]
                
                closest_term.sort(key=lambda ct: ct[1])
                
                if not len(closest_term):
                    
                    continue
                
                else:
                    
                    closest_term, closest_term_distance = closest_term[0]
                    
                    # If closest term is within max_distance from cluster_parent, then add to cluster
                    # Else, create new leaf
                    
                    if closest_term_distance <= self.max_distance:
                    
                        clusters[closest_term] = l_clust
                        
                        clusters_parents[l_clust] = closest_term
                        
                        break
                    
                    else:
                        
                        leaves.append(closest_term)
                        
                        new_cluster_n = max(clusters.values()) + 1
                        
                        clusters[closest_term] = new_cluster_n
                        
                        clusters_parents[new_cluster_n] = closest_term
                        
                        break
            
            if clusters == old_clusters:
                
                toggle = False
            
        log.append(f'Created {len(leaves) - original_leaves_num} new leaves')
            
        log.append(f'Found {max(clusters.values()) + 1} clusters')
        
        # Keep shared terms that don't have any children among other shared terms
        
        shared_terms = to_be_removed
        
        shared_terms = [s for s in shared_terms if np.isin(list(self.get_children_terms([s])), shared_terms).sum() == 1]
        
        # Simplify clusters using shared parent terms previously removed
        
        for s in shared_terms:
            
            # Find clusters sharing s
            
            s_children = self.get_children_terms([s])
        
            clusters_to_merge = [cl for cl,h in clusters_parents.items() if h in s_children]
            
            # Only keep those whose distance from s is <= max_distance
            
            clusters_to_merge = [cl for cl in clusters_to_merge if self.get_parent_child_distance(s, clusters_parents[cl]) <= self.max_distance]
            
            if len(clusters_to_merge):
                
                new_cluster_n = min(clusters_to_merge)
                
                # Update clusters
                
                elements_to_update = [h for h,cl in clusters.items() if cl in clusters_to_merge] + [s]
                
                for e in elements_to_update:
                    
                    clusters[e] = new_cluster_n
                
                # Update clusters_parents
                
                for cl in clusters_to_merge:
                    
                    if cl == new_cluster_n:
                        
                        clusters_parents[new_cluster_n] = s
                    
                    else:
                        
                        _ = clusters_parents.pop(cl, None)
                
                # Log
                
                log.append(f'Merged clusters {clusters_to_merge} into cluster {new_cluster_n}')
        
        # Add unclustered terms as individual clusters
        
        for term,cl in clusters.items():
            
            if cl == -1:
                
                new_cl = max(clusters.values()) + 1
                
                clusters[term] = new_cl
                
                clusters_parents[new_cl] = term
        
        # Re-index clusters
        
        reindexing = {old : new for new,old in enumerate(np.sort(list(clusters_parents.keys())))}
        
        clusters = {h : reindexing[cl] for h,cl in clusters.items() if cl != -1}
        
        clusters_parents = {reindexing[cl] : h for cl,h in clusters_parents.items()}
    
        log.append(f'Reduced data to {max(clusters.values()) + 1} clusters')
        
        self.clusters = clusters
        
        self.clusters_parents = clusters_parents
        
        self.log = log
    
    ### ------------------------------------ ###
    ### UTILS                                ###
    ### ------------------------------------ ###
    
    def get_children_terms(self, c):
        
        children = set(self.parent_child.loc[self.parent_child.parent.isin(c), 'child'].to_list())
        
        children.update(c)
    
        if children == c:
            
            return c
        
        else:
            
            return self.get_children_terms(children)
    
    ### ------------------------------------ ###
    
    def get_parent_terms(self, p):
        
        parents = set(self.parent_child.loc[self.parent_child.child.isin(p), 'parent'].to_list())
        
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
            
            children = set(self.parent_child.loc[self.parent_child.parent.isin(children), 'child'].to_list())
        
        return n
    
    ### ------------------------------------ ###
    
    def find_closest_common_ancestor(self, element_1, element_2):
        
        # Get ancestors in common
        
        parents_1 = self.get_parent_terms([element_1])
        
        parents_2 = self.get_parent_terms([element_2])
        
        common_ancestors = [p1 for p1 in parents_1 if p1 in parents_2 and p1 != '']
        
        if not len(common_ancestors):
            
            return '', -1
        
        else:
            
            # Find shortest distance
        
            elements_distance = []
            
            for ca in common_ancestors:
                
                distance_1 = self.get_parent_child_distance(ca, element_1)
                
                distance_2 = self.get_parent_child_distance(ca, element_2)
                
                elements_distance.append((ca, distance_1 + distance_2 + 1))
            
            elements_distance.sort(key=lambda ed: ed[1])
            
            # Closest ancestor
            
            closest_ancestor, branch_length = elements_distance[0]
            
            return closest_ancestor, branch_length
    
    ### ------------------------------------ ###
    
    def save_data(self):
        
        with open(f'{self.output_prefix}_simplified.txt', 'w') as term_out:
            
            term_out.write('\n'.join(self.simplified_terms))

        with open(f'{self.output_prefix}_clusters.txt', 'w') as clusters_out:
            
            clusters_out.write('\n'.join(['term\tcluster'] + [f'{h}\t{cl}' for h,cl in self.clusters.items()]))

        with open(f'{self.output_prefix}_clustering_log.txt', 'w') as log_out:
            
            log_out.write('\n'.join(self.log))

### ---------------------------------------- ###

import numpy as np
import pandas as pd
