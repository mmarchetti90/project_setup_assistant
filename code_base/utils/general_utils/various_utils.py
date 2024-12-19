#!/usr/bin/env python3

"""
Collection of useful functions
"""

### ---------------------------------------- ###
### UTILS                                    ###
### ---------------------------------------- ###

def get_points_in_between(start, stop):

    """
    Given two points, it finds the necessary ones to connect them.
    """
    
    distance = round(((stop[0] - start[0])**2 + (stop[1] - start[1])**2)**0.5)
    x_shift, y_shift = stop - start
    
    points = []
    
    for d in range(distance):
        
        new_point = [start[0] + round(d * x_shift / distance), start[1] + round(d * y_shift / distance)]
        points.append(new_point)
    
    if list(stop) not in points:
        
        points.append(stop)
    
    points = np.array(points)
    
    return points

### ---------------------------------------- ###

def kneedle(vector, sort_vector=True):
    
    """
    Kneedle to find threshold cutoff.
    """
    
    if sort_vector:
        
        vector = np.sort(vector)[::-1]
    
    # Find gradient and intercept
    x0, x1 = 0, len(vector)
    y0, y1 = max(vector), min(vector)
    gradient = (y1 - y0) / (x1 - x0)
    intercept = y0
    
    # Compute difference vector
    difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(vector)]
    
    # Find max of difference_vector and define cutoff
    cutoff_index = difference_vector.index(max(difference_vector))
    cutoff_value = vector[cutoff_index]
    
    return cutoff_index, cutoff_value

### ---------------------------------------- ###

def running_average(vector, window=10, valid_mode=False):

    """
    Computes the running average of a list of values or numpy array given a window.
    """
    
    half_window = window // 2
    indexes = np.arange(half_window, len(vector) - half_window, 1)
    smooth = np.array([np.mean(vector[i - half_window : i + half_window]) for i in indexes])
    
    if not valid_mode:
        
        tail_left_indexes = np.arange(0, indexes[0], 1)
        tail_left_smooth = np.array([np.mean(vector[max(0, i - half_window) : i + half_window]) for i in tail_left_indexes])
        tail_right_indexes = np.arange(indexes[-1] + 1, len(vector), 1)
        tail_right_smooth = np.array([np.mean(vector[i - half_window : min(len(vector), i + half_window)]) for i in tail_right_indexes])
        
        indexes = np.concatenate([tail_left_indexes, indexes, tail_right_indexes])
        smooth = np.concatenate([tail_left_smooth, smooth, tail_right_smooth])
    
    return indexes, smooth

### ---------------------------------------- ###

def sort_coords(pnts):
    
    """
    Sort coordinates that make up a line.
    """
    
    # Run KNN
    nbrs = NearestNeighbors(n_neighbors=3, algorithm='kd_tree').fit(pnts)
    indices = nbrs.kneighbors(pnts, return_distance=False)
    
    # Init sorted list of indices
    sorted_idx = indices[0, [1, 0, 2]].tolist()
    
    # Add points right
    toggle = True
    while toggle:
        
        i = sorted_idx[-1]
        new_idx = indices[i, 1:]
        if new_idx[0] not in sorted_idx:
            
            next_i = new_idx[0]
        
        elif new_idx[1] not in sorted_idx:
            
            next_i = new_idx[1]
            
        else:
            
            toggle = False
            continue
        
        sorted_idx.append(next_i)
    
    # Invert sorted_idx
    sorted_idx = sorted_idx[::-1]
    
    # Add points left
    toggle = True
    while toggle:
        
        i = sorted_idx[-1]
        new_idx = indices[i, 1:]
        if new_idx[0] not in sorted_idx:
            
            next_i = new_idx[0]
        
        elif new_idx[1] not in sorted_idx:
            
            next_i = new_idx[1]
            
        else:
            
            toggle = False
            continue
        
        sorted_idx.append(next_i)
    
    sorted_pnts = pnts[sorted_idx]
    
    return sorted_pnts

### ------------------MAIN------------------ ###

try:
    
    import numpy as np
    
    from sklearn.neighbors import NearestNeighbors
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()