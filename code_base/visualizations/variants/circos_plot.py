#!/usr/bin/env python3

"""
Example of circos plot of chromosomal aberrations by patient

INPUT
-----

patient_data
    Tab-separated table with shape N * (M + 1), with N being patients and M the chromosomal
    aberrations, one-hot encoded. THe first column is unique patients identifiers.
    
    e.g.
              Normal  del5q  del7q  del9q
    Patient1  1       0      0      0
    Patient2  0       1      0      1
    ...
    PatientN  0       0      1      0
    
"""

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

import pandas as pd

from pycirclize import Circos
from sys import argv

### Load data

data_path = argv[argv.index('--patient_data') + 1]

data = pd.read_csv(data_path, sep='\t', index_col=0)

### Sort data

codings = [(n, ''.join(row.astype(str))) for n,row in data.iterrows()]
codings.sort(key=lambda c: c[1], reverse=True)
data = data.loc[[c[0] for c in codings],]

### Circos plot

# Init circos object
circos = Circos({'' : data.shape[0]}, start=90, end=360)

# Define sector properties
sector = circos.sectors[0]
sector.axis(fc='none', ls='-', lw=2, ec='black', alpha=1)

# Add tracks
bottom_gap = 10
track_size = (100 - bottom_gap) / (data.shape[1] + 1)
text_common_kws = dict(ha='left', va='center', size=16)
#colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:cyan', 'tab:pink', 'tab:gray'] # Named colors for 6 samples
for n,col in enumerate(data.columns.values):
    
    track = sector.add_track((100 - track_size * (n + 1), 100 - track_size * n))
    track.axis(fc='none', ec='black', lw=0.5)
    #track.bar(range(data.shape[0]), data[col], width=1, color=colors[n], align='edge', lw=0.25)
    track.bar(range(data.shape[0]), data[col], width=1, color='red', align='edge', lw=0.25)
    
    circos.text(f' {col}', r=(100 - track_size * (n + 0.5)), color='black', **text_common_kws)

# Plot
#circos.savefig('chromosomal_aberrations_by_patient.png', dpi=300)
circos.savefig('chromosomal_aberrations_by_patient.svg', dpi=300)
#fig = circos.plotfig() # For testing
