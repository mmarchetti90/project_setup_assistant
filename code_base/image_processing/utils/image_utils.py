#!/usr/bin/env python3

"""
Collection of useful functions for image processing
"""

### ---------------------------------------- ###
### IMAGE THRESHOLDING                       ###
### ---------------------------------------- ###

def threshold_movie(movie, output_prefix='sample'):

    """
    Measures intensities along the 4 axis of the image, then pools the measurements, sorts them, and
    uses a Kneedle approach to find an optimal threshold.
    """

    print('Thresholding')

    # Find thresholds via kneedle
    thresholds = []
    for t,m in enumerate(movie):
        
        threshold = threshold_frame(m)
        
        thresholds.append(threshold)
    
    # Running average
    _, thresholds = running_average(thresholds, int(len(thresholds) * 0.2), valid_mode=False)

    # Save thresholds to file and plot data
    save_threshold_diagnostics(thresholds, f'{output_prefix}_thresholds')
    
    return thresholds

### ---------------------------------------- ###

def threshold_frame(frame):

    """
    Measures intensities along the 4 axis of the image, then pools the measurements, sorts them, and
    uses a Kneedle approach to find an optimal threshold.
    """

    print('Thresholding')

    # Get profile along several axes
    ax_1 = frame[frame.shape[0] // 2,]
    ax_2 = frame[:, frame.shape[1] // 2]
    ax_3_pnts = get_points_in_between(np.array([0, 0]), np.array(frame.shape) - 1)
    ax_3 = frame[ax_3_pnts[:,0], ax_3_pnts[:,1]]
    ax_4_pnts = get_points_in_between(np.array([frame.shape[0] - 1, 0]), np.array([0, frame.shape[1]]) - 1)
    ax_4 = frame[ax_4_pnts[:,0], ax_4_pnts[:,1]]
    
    # Concatenate
    ax = np.concatenate([ax_1, ax_2, ax_3, ax_4])
    
    # Smoothing
    averaging_window = 10
    _, ax = running_average(ax, averaging_window)
    
    # Sorting
    ax = np.sort(ax)[::-1]
    
    # Kneedle on ax
    _, threshold = kneedle(ax)
    
    return threshold

### ---------------------------------------- ###

def save_threshold_diagnostics(thrs, prefix='threshold'):

    # Save to file
    data = pd.DataFrame({'frame' : range(len(thrs)),
                         'threshold' : thrs})
    data.to_csv(f'{prefix}.tsv', sep='\t', index=False, header=True)
    
    # Smooth
    idx, smooth_thrs = running_average(thrs, 10)
    
    # Plot data
    plt.figure(figsize=(15, 3))
    plt.plot(thrs, 'g', lw=3)
    plt.plot(idx, smooth_thrs, 'r', lw=3)
    plt.xlabel('Time (frames)')
    plt.tight_layout()
    plt.savefig(f'{prefix}.png', dpi=300)
    plt.close()

### ---------------------------------------- ###
### MOVIE MASKING                            ###
### ---------------------------------------- ###

def mask_movie(movie, thresholds, size_filter=True, min_size=100, output_prefix='sample'):

    print('Masking')

    masks = []
    for t,(m,threshold) in enumerate(zip(movie, thresholds)):
    
        # Initial mask
        mask = (m > threshold)
    
        # Clean mask by erosion + expansion (this removes most small objects)
        mask = maxfil(minfil(mask, 3), 3)
    
        # Label foreground objects
        labels, labels_num = label(mask, np.ones((3, 3)))
        labels_areas = [(lb + 1, (labels == lb + 1).sum()) for lb in range(labels_num)]
        labels_areas.sort(key=lambda x: x[1], reverse=True)
    
        # Filter ojects by size
        if labels_num > 5 and size_filter:
    
            index, _ = kneedle([la[1] for la in labels_areas])
    
        else:
        
            index = labels_num
        
        good_labels = [la[0] for la in labels_areas[:index] if la[1] > min_size]
    
        # Redefine mask
        mask = np.ones(m.shape)
        mask[np.where(np.isin(labels, good_labels))] = 0
        
        # Fill holes
        for _ in range(3):
            
            mask = maxfil(minfil(mask, 5), 5)
        
        masks.append(mask)
        
        print(f'\rMasked {t + 1} / {movie.shape[0]} frames', end='')
    
    # Make diagnostic mask for m
    diagnostics = []
    for mk,m in zip(masks, movie):
        
        outlines = 255 * ((mk != 0).astype(int) - minfil(mk != 0, 3).astype(int))
        diag = np.stack([outlines, m, np.zeros(m.shape, dtype=int)], axis=-1)
        diagnostics.append(diag)
    
    # Converting to array
    masks = np.array(masks)
    diagnostics = np.array(diagnostics)

    # Save mask to file
    save_tiff(masks, f'{output_prefix}_mask')
    save_tiff(diagnostics, f'{output_prefix}_mask-diagnostics')

    mask = masks.copy()
    
    return mask

### ---------------------------------------- ###
### MOVIE BACKGROUND REMOVAL                 ###
### ---------------------------------------- ###

def clean_movie(mov, background_multiplier=0.5):
    
    """
    Cleans the movie using a no-neigbor deblurring algorithm.
    """

    # Using no-neighbour deblurring (i.e. simulating background by blurring image using different sigmas) to find foreground
    print("Running no-neighbour deblurring")
    
    # Smoothening the image using a mean filter
    clean_mov = correlate(mov, np.ones((1, 3, 3)) / 9, mode="same")
    
    # To make the calculation faster, mov is downscaled to find background
    clean_mov = np.array([np.array(img.fromarray(t).resize((int(t.shape[1] / 4), int(t.shape[0] / 4)))) for t in mov])
    
    # Simulating background
    background = np.zeros(clean_mov.shape)
    
    for sigma in [5, 10, 20]:
        
        background = np.max((background, gaussian(clean_mov, (0, sigma, sigma), mode="reflect")), 0)
    
    clean_mov = (clean_mov - background * background_multiplier) # Adjusting background intensity

    # Reshaping to 8bit
    clean_mov = (clean_mov * 255 / np.max(clean_mov)).astype('int16')
    
    # Restore original size
    clean_mov = np.array([np.array(img.fromarray(m).resize((mov.shape[2], mov.shape[1]))) for m in clean_mov])
    
    clean_mov[clean_mov < 0] = 0

    return clean_mov

### ------------------------------------ ###
### IMAGE IMPORT/EXPORT                  ###
### ------------------------------------ ###

def image_importer(path, max_frames=10000):

    """
    Imports an image as a Numpy array.
    """
    
    print(f'Importing image {path}')
    
    raw = img.open(path)
    
    movie_arr = []
    
    for r in range(0, max_frames):
        
        try:
            
            raw.seek(r)
            movie_arr.append(np.array(raw))
        
        except:
            
            break

    movie_arr = np.array(movie_arr, dtype='int16') # int8 was giving problems, i.e. bright pixels, close to 256, were converted to 0
    
    raw.close()
    
    print("Image has been successfully imported")
    
    return movie_arr

### ------------------------------------ ###

def save_tiff(img_arr, output_prefix='pic'):
    
    """
    Exporting mask as multipage tiff
    """
    
    if img_arr[0].shape[-1] == 3: # RGB image
    
        img_arr = [img.fromarray(ia.astype(np.uint8), 'RGB') for ia in img_arr]
    
    else: # Mono-channel image
        
        img_arr = [img.fromarray(ia.astype(np.uint8)) for ia in img_arr]
    
    img_arr[0].save(f'{output_prefix}.tif', "TIFF", save_all=True, append_images=img_arr[1:])

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
    import pandas as pd
    
    from matplotlib import pyplot as plt
    from PIL import Image as img
    #from scipy.ndimage import binary_fill_holes as fill
    from scipy.ndimage import gaussian_filter as gaussian
    from scipy.ndimage import label
    from scipy.ndimage import minimum_filter as minfil
    from scipy.ndimage import maximum_filter as maxfil
    from scipy.signal import correlate #Easier math than convolution
    from sklearn.neighbors import NearestNeighbors
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()