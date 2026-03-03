#!/usr/bin/env python3

"""
Based on "Analyzing nested experimental designs-A user-friendly resampling method to determine experimental significance"
Kulkarni et al
PLoS Comput Biol, 2022
https://pubmed.ncbi.nlm.nih.gov/35500032/
"""

### ---------------------------------------- ###

class nested_analysis:
    
    """
    Class for the analysis of Nested Experimental Designs.
    Based on Kulkarni et al, PLOS Computational Biology, 2022
    See https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010061

    FIX THE FOLLOWING DESCRIPTION
    
    Parameters
    ----------
    data : Pandas data frame
        Data to be processed.
    treatment_col : string
        Column in data with treatment information.
    biological_replicate_col : string
        Column in data with biological replicates information.
    measurement : string
        Column in data with the measurement from biological and technical replicates.
    bootstrapping_iterations : int, optional
        Number of bootstrapping iterations.
        Default=20
    max_permutations : int, optional
        Maximum number of permutations on bootstrapped data.
        Default=1000

    Methods
    -------
    analyze()
        Analyzes data and returns a two-tailed p-value.
    """
    
    def __init__(self, data, treatment_col, biological_replicate_col, measurement, bootstrapping_iterations=20, max_permutations=1000):
        
        self.data = data
        self.tr_col = treatment_col
        self.tr_class = list(set(data.loc[:, self.tr_col].values))
        self.bio = biological_replicate_col
        self.m = measurement
        self.boots = bootstrapping_iterations
        
        # Get samples and sample size to compute the maximum permutations
        self.samples = list(set(self.data[self.bio].values))
        self.cl1_size = np.unique(self.data.loc[self.data[self.tr_col] == self.tr_class[0], self.bio]).shape[0]
        self.cl2_size = np.unique(self.data.loc[self.data[self.tr_col] == self.tr_class[1], self.bio]).shape[0]
        max_possible_perms = self.factorial(len(self.samples)) / (self.factorial(self.cl1_size) * self.factorial(len(self.samples) - self.cl1_size))
        self.perms = min(max_possible_perms, max_permutations)
    
    ### ------------------------------------ ###
    ### STATISTIC TESTING                    ###
    ### ------------------------------------ ###
    
    def analyze(self):
        
        # Compute statistic for data
        print('Computing statistics for base data')
        cl1_dat = self.data.loc[self.data[self.tr_col] == self.tr_class[0],].groupby(by=self.bio, axis=0).mean(numeric_only=True)[self.m].values
        cl2_dat = self.data.loc[self.data[self.tr_col] == self.tr_class[1],].groupby(by=self.bio, axis=0).mean(numeric_only=True)[self.m].values
        t, _ = ttest_ind(cl1_dat, cl2_dat, equal_var=False, alternative='two-sided')
        self.base_t = t
        
        # Bootstrapping samples and null distribution of t values
        null_distribution = []
        for _ in range(self.boots):
            
            # Bootstrapping
            bootstrapped_means = {s : self.bootstrap_sample(self.data.loc[self.data[self.bio] == s, self.m]) for s in self.samples}
            
            # Permuting
            for n,comb in enumerate(combinations(self.samples, self.cl1_size)):
                
                if n >= self.perms:
                    
                    break
                
                cl1_keys = comb
                cl2_keys = [s for s in self.samples if s not in cl1_keys]
                t, _ = ttest_ind([bootstrapped_means[s] for s in cl1_keys],
                                 [bootstrapped_means[s] for s in cl2_keys],
                                 equal_var=False, alternative='two-sided')
                
                null_distribution.append(t)
        
        self.null_distribution = np.array(null_distribution)
        
        # Final statistics
        p = (((self.null_distribution < - abs(self.base_t)).sum() +
              (self.null_distribution > abs(self.base_t)).sum())
             / len(self.null_distribution))
        self.nestest_pval = p
        
        # Diagnostic plot
        plt.figure(figsize=(10, 5))
        (y, x, _) = plt.hist(null_distribution, bins=100)
        plt.vlines(self.base_t, 0, y.max(), 'black', '--')
        plt.vlines(- self.base_t, 0, y.max(), 'black', '--')
        plt.xlabel('t-distribution')
        plt.ylabel('Count')
    
    ### ------------------------------------ ###
    ### UTILS                                ###
    ### ------------------------------------ ###
    
    def factorial(self, val):
        
        if val > 1:
            
            return (val * self.factorial(val - 1))
        
        else:
            
            return 1
    
    ### ------------------------------------ ###
    
    @staticmethod
    def bootstrap_sample(values):
        
        n = len(values)
        bootstrapped = np.random.choice(values, n, replace=True)
        
        return np.mean(bootstrapped)

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import seaborn as sns

from itertools import combinations
from matplotlib import pyplot as plt
from scipy.stats import ttest_ind
