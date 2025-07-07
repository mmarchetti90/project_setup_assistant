#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import numpy as np

### CLASSES AND FUNCTIONS ------------------ ###

def split_data(y: np.ndarray, n_splits: int=10, label_batch_size: int=5, shuffle: bool=False) -> (list, list):
    
    """
    KFold-like split for unbalanced labels.
    Will make sure that the same number of samples are used for all labels and that all data is
    used at least once.
    Works by STUB
    
    Parameters
    ----------
    y : np.array()
        Numpy array of labels
    n_splits : int, optional
        Number of dataset splits
    label_batch_size : int, optional
        Number of samples from each label to be included in each split
    shuffle : bool, optional
        Set to True if labels should be shuffled prior to split
    """
    
    # Sanity check
    
    y = np.array(y).ravel()
    
    # Find all unique labels
    
    unique_labels = np.sort(np.unique(y))
    
    labels_count = [(l, (y == l).sum()) for l in unique_labels]
    
    smallest_label_count = min(labels_count, key=lambda lc: lc[1])
    
    largest_label_count = max(labels_count, key=lambda lc: lc[1])
    
    if smallest_label_count[1] <= label_batch_size:
        
        print(f'ERROR: smallest label ({smallest_label_count[0]}) has too few samples ({smallest_label_count[1]}) for label_batch_size={label_batch_size}')
        
        return [], []
    
    # Labels indexes
    
    idx = {l : np.where(y == l)[0] for l in unique_labels}
    
    # Shuffle
    
    if shuffle:
        
        idx = {l : np.random.choice(i, size=i.shape[0], replace=False) for l,i in idx.items()}
    
    # Splits
    
    min_splits = (largest_label_count[1] // label_batch_size) + (1 if (largest_label_count[1] % label_batch_size) > 0 else 0)
    
    n_splits = n_splits if n_splits >= min_splits else min_splits
    
    splits = {s : [] for s in range(n_splits)}
    
    for l,i in idx.items():
        
        # Extend indexes array to reach a size of n_splits * label_batch_size

        original_i_size = i.shape[0]

        while i.shape[0] < n_splits * label_batch_size:
            
            i = np.concatenate([i, np.random.choice(i[:original_i_size], size=original_i_size, replace=False)])
        
        for s in range(n_splits):
            
            start, end = int(s * label_batch_size), int((s + 1) * label_batch_size)
            
            splits[s].extend(i[start : end].tolist())
        
    # Divide into train and validation
    
    train_splits = list(splits.values())
    
    validation_splits = []
    
    for snum in range(n_splits):
        
        # Choose random split to serve as validation set
        
        val_s = np.random.choice([s for s in range(n_splits) if s not in validation_splits + [snum]], size=1)[0]
        
        validation_splits.append(val_s)
        
    validation_splits = [splits[vs] for vs in validation_splits]
    
    return train_splits, validation_splits

### ---------------------------------------- ###

if __name__ == '__main__':
    
    pass
