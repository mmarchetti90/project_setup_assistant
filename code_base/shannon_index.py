#!/usr/bin/env python3

### ---------------------------------------- ###

def compute_shannon_index(species_abundance):

    """
    Function for computing Shannon's diversity index

    Parameters
    ----------
    species_abundance : list or np.array
        Abundance of the individual species in the sample
    """

    species_abundance = np.array(species_abundance)

    tot = species_abundance.sum()
    prop = species_abundance["abundance"] / tot
    ln_prop = np.log(prop)
    s_index = - sum(prop * ln_prop)

    return s_index

### ------------------MAIN------------------ ###

import numpy as np
