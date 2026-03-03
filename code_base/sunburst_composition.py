#!/usr/bin/env python3

"""
Plotting bacteria composition in replicates of sample groups for genus and species
Sunbursts are organized in a grid with shape len(sample groups) * len(samples)

INPUTS
------

data
    Tab-delimited table with the following necessary columns:
        sample_id : unique sample identifiers
        group_id : unique sample group identifiers
        genus : genus information
        species : species information. Set to "sp." if unknown
        abundance : species abundance in each sample
"""

### ---------------------------------------- ###

def rgb2hex(rgb):
    
    r, g, b = rgb
    r, g, b = int(255 * r), int(255 * g), int(255 * b)
    
    return f"#{r:02x}{g:02x}{b:02x}"

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from matplotlib import pyplot as plt
from plotly.subplots import make_subplots
from sys import argv

### Load data

data_file = argv[argv.index('--data') + 1]

data = pd.read_csv(data_file, sep='\t')

data.fillna(0, inplace=True)

### Genus and species level sunburst

all_genus = np.sort(np.unique(data.genus))
color_palette = ([rgb2hex(plt.get_cmap('tab20')(n)[:3]) for n in range(20)] +
                 [rgb2hex(plt.get_cmap('tab20')(n)[:3]) for n in range(20)] +
                 [rgb2hex(plt.get_cmap('tab20')(n)[:3]) for n in range(20)] +
                 [rgb2hex(plt.get_cmap('tab20')(n)[:3]) for n in range(20)])
genus_color = {g : color_palette[n * 2] for n,g in enumerate(all_genus)}
species_color = {g : color_palette[n * 2 + 1] for n,g in enumerate(all_genus)}

groups = np.unique(data.group.values)
n_groups = len(groups)

samples = {g : np.unique(data.loc[data.group == g, 'sample_id'].values) for g in groups}
n_samples = max([len(s) for s in samples.values()])

#subplot_titles = [f'<b>{s}</b>' for g_samples in samples.values() for s in g_samples]
subplot_titles = [f'<b>rep_{ns}</b>' for ns in range(n_samples)] + ['' for _ in range(n_samples * (n_groups - 1))]
row_labels = [f'<b>{g}</b>' for g in groups]

fig = make_subplots(rows=n_groups,
                    cols=n_samples,
                    specs=[[{'type': 'sunburst'} for _ in range(n_samples)] for _ in range(n_groups)],
                    subplot_titles=subplot_titles,
                    row_titles=row_labels,
                    horizontal_spacing=0.01,
                    vertical_spacing=0.01)

for ng,gr in enumerate(groups):
    
    g_samples = samples[gr]
    
    for ns,s in enumerate(g_samples):

        sub_data = data.loc[(data.group == gr) & (data.sample_id == s), ['genus', 'species', 'abundance']]
        sub_data.replace({'species' : {'sp.' : ''}}, inplace=True)
        
        #sub_data['abundance'] = np.round(sub_data['abundance'])
        
        genus = np.sort(np.unique(sub_data.genus))
        color_map_sub = {g : c for g,c in genus_color.items() if g in genus}
        
        # Add root
        labels, parents, ids, values, colors = [], [], [], [], []
        
        for g in genus:
            
            g_sub = sub_data.loc[sub_data.genus == g,]
        
            # Add genus info
            labels.append(g.replace(' ', '<br>'))
            parents.append('')
            ids.append(g)
            values.append(int(g_sub['abundance'].sum()))
            colors.append(genus_color[g])
            
            for _,(_,sp,v) in g_sub.iterrows():
                
                if len(sp):
                    
                    labels.append(sp.replace(' ', '<br>'))
                    parents.append(g)
                    ids.append(f'{g} - {sp}')
                    values.append(int(v))
                    colors.append(species_color[g])
        
        fig.add_trace(go.Sunburst(ids=ids,
                                  labels=labels,
                                  parents=parents,
                                  values=values,
                                  branchvalues='total',
                                  insidetextfont=dict(color='#000000', size=30),
                                  insidetextorientation='radial',
                                  marker=dict(line=dict(width=2, color='#000000'),
                                              colors=colors),
                                  domain=dict(x=[(1 / n_samples) * ns,
                                                 (1 / n_samples) * (ns + 1)],
                                              y=[(1 / n_groups) * ng,
                                                 (1 / n_groups) * (ng + 1)])),
                      row=ng + 1, col=ns + 1)

fig.update_layout(margin=dict(t=150, l=200, r=0, b=0))
fig.for_each_annotation(lambda a: a.update(x=-0.015, textangle=-90) if a.text in row_labels
                        else a.update(y=1.005) if a.text in subplot_titles
                        else())
fig.update_annotations(font_size=100, font=dict(color='black')) # Subplot titles size

#fig.update_layout(margin=dict(t=0, l=0, r=0, b=0),
#                  uniformtext=dict(minsize=20, mode='hide'))

fig.write_image('genus_composition_sunburst.pdf', width=1048 * n_samples, height=1048 * n_groups)
