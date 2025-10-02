#!/usr/bin/env python3

"""
Class for finding outlier values between two samples

Parameters
----------
data : pandas frame
    Pandas dataframe with 3 columns: variable identifier, sample 1 values, and sample 2 values
normalize_data : bool
    Set to True to normalize data
"""

### ---------------------------------------- ###

class find_outliers:

    def __init__(self, data, normalize_data=True):
        
        data[data.columns.values[1:]] = data.loc[:, data.columns.values[1:]].astype(float)
        
        if normalize_data:
            
            data = self.normalize_counts(data.copy())
        
        self.variable_names = data.iloc[:, 0].values.ravel()
        self.x_name, self.y_name = data.columns[1:]
        self.x_val, self.y_val = data.iloc[:, 1].values.ravel(), data.iloc[:, 2].values.ravel()
        
    ### ------------------------------------ ###
    ### NORMALIZE DATA                       ###
    ### ------------------------------------ ###

    @staticmethod
    def normalize_counts(raw_counts):
        
        normalized_counts = raw_counts.copy()
            
        # Proportional fitting
        library_sizes = normalized_counts.iloc[:, 1:].sum(axis=0).values
        median_size = np.median(library_sizes)
        normalized_counts.iloc[:, 1:] = normalized_counts.iloc[:, 1:].div(library_sizes / median_size, axis=1)
        
        # log1p
        normalized_counts.iloc[:, 1:] = np.log1p(normalized_counts.iloc[:, 1:])
        
        # Proportional fitting
        library_sizes = normalized_counts.iloc[:, 1:].sum(axis=0).values
        median_size = np.median(library_sizes)
        normalized_counts.iloc[:, 1:] = normalized_counts.iloc[:, 1:].div(library_sizes / median_size, axis=1)
        
        return normalized_counts
    
    ### ------------------------------------ ###
    ### OUTLIER DETECTION                    ###
    ### ------------------------------------ ###
    
    def get_outliers(self, n_bins=20, std_max_distance=2):
        
        """
        Function to find outliers based on the distance from a trendline
        
        Steps
        -----
        1. Correlate x and y values to define a trendline
        2. For each point, determine the closest one on the trendline
        3. Bin values based on x range
        4. For each bin, find a set of point that are the closest to the trendline section, then correlate the values of the two samples to define a local trendline
        3. Mark samples as outliers if their distance from the trendline is > 2 * standard deviation
        
        Parameters
        ----------
        n_bins : int
            Number of bins for local clustering of data points
        std_max_distance : int
            Number of standard deviations of distances of points from a trendline
            Used to define outliers
        """
        
        # Interpolate x and y to find a gross trendline, then for each x,y pair, find closest point on trendline
    
        slope, intercept, r, p, closest_x, closest_y = self.find_closest_point_on_trendline(self.x_val, self.y_val)
        
        y_predict = slope * self.x_val + intercept
        
        # Bin along the x axis
        
        xmin, xmax = self.x_val.min(), self.x_val.max()
        bin_size = (xmax - xmin) / (n_bins - 1)
        
        # For each bin, find outliers
        
        outliers_idx = []
        
        for n in range(n_bins):
            
            # Subset values
            
            x0, x1 = xmin + bin_size * n, xmin + bin_size * (n + 1)
            
            subset_idx = np.where((closest_x >= x0) & (closest_x < x1))[0]
            
            x_val_sub, y_val_sub = self.x_val[subset_idx], self.y_val[subset_idx]
            
            closest_x_sub, closest_y_sub = closest_x[subset_idx], closest_y[subset_idx]
    
            # Correlate x and y
            
            try:
                
                # Find points distance to trendline
                
                distances_to_trendline = [np.sqrt((xa - xb)**2 + (ya - yb)**2) for xa,ya,xb,yb in zip(x_val_sub, y_val_sub, closest_x_sub, closest_y_sub)]
                
                std = np.std(distances_to_trendline)
        
                # Define outliers based on > 2 std from best-fit line
                
                outliers_idx += subset_idx[np.where(distances_to_trendline > (std_max_distance * std))[0]].tolist()
            
            except:
                
                pass
    
        # Diagnostic plot
        
        plt.figure(figsize=(10, 10))
        plt.xlabel(f'{self.x_name} counts')
        plt.ylabel(f'{self.y_name} counts')
        plt.scatter(self.x_val, self.y_val, color='lightgray', lw=0.5, edgecolor='black')
        plt.scatter(self.x_val[outliers_idx], self.y_val[outliers_idx], color='red', lw=0.5, edgecolor='black')
        plt.plot(self.x_val, y_predict, 'green')
        plt.savefig('outliers.png', dpi=300)
        plt.close()
    
        self.outliers = np.array(self.variable_names)[outliers_idx]
    
    ### ------------------------------------ ###
    
    @staticmethod
    def find_closest_point_on_trendline(x_val, y_val):
        
        # Interpolate x and y to find a trendline
    
        slope, intercept, r, p, std_err = linregress(x_val, y_val, alternative='two-sided')
        
        # For each x,y pair, find closest point on trendline
        
        closest_x, closest_y = [], []
        
        slope_perpendicular = - 1 / slope
        
        for x,y in zip(x_val, y_val):
            
            intercept_perpendicular = y - slope_perpendicular * x
            
            x_intercept = (intercept_perpendicular - intercept) / (slope - slope_perpendicular)
        
            y_intercept = slope * x_intercept + intercept
            
            #y_closest = slope_perpendicular * x_closest + intercept_perpendicular
            
            closest_x.append(x_intercept)
            
            closest_y.append(y_intercept)
        
        closest_x, closest_y = np.array(closest_x), np.array(closest_y)
        
        return slope, intercept, r, p, closest_x, closest_y

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from scipy.stats import linregress

# Load counts data

counts_path = '../../../2_processed_data/rna_data/2_gene_counts/MergedGeneCounts_Firststrand.tsv'
counts = pd.read_csv(counts_path, sep='\t')
counts[counts.columns.values[1:]] = counts.loc[:, counts.columns.values[1:]].astype(float)

# Remove undetected in both samples

counts = counts.loc[counts.iloc[:, 1:].sum(axis=1) > 0,]










