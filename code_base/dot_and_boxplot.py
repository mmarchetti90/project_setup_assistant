#!/usr/bin/env python3

### ---------------------------------------- ###

def dot_and_boxplot(data, x, y, hue=None, dodge=False, xlab='X', ylab='Y', size=10, linewidth=1, edgecolor='black', palette='tab10', figsize=(10, 5), plot_type='swarm', figname='figure.png'):

    """
    Function for plotting dots and an overlaid boxplot

    Parameters
    ----------
    data : Pandas DataFrame
        Dataframe to plot
    x : string
        Column for x-axis
    y : string
        Column for y-axis
    hue : string, optional
        Column for samples' color
        Default=None
    dodge : bool, optional
        If True, separates samples based on hue
        Default=False
    xlab : string, optional
        Label for the x-axis
        Default='X'
    ylab : string, optional
        Label for the y-axis
        Default='Y'
    size : int, optional
        Size of dots
        Default=10
    linewidth : int, optional
        Size of dots edge
        Default=1
    edgecolor : string, optional
        Color of dots edge
        Default='black'
    palette : string or list of colors, optional
        Palette for dots
        Default='tab10'
    figsize : tuple of ints, optional
        Size of output figure
        Default=(10, 5)
    plot_type : string, optional
        Use 'swarm' for swarmplot and 'strip' for stripplot
    figname : string, optional
        Name of output file
        Default='figure.png'
    """

    plt.figure(figsize=figsize)

    if plot_type == 'swarm':

        dots = sns.swarmplot(data_grouped, x=x, y=y, hue=hue, dodge=dodge,
                             size=size, linewidth=linewidth, edgecolor=edgecolor, palette=colors)

    elif plot_type == 'swarm':

        dots = sns.stripplot(data_grouped, x=x, y=y, hue=hue, dodge=dodge,
                             size=size, linewidth=linewidth, edgecolor=edgecolor, palette=colors)

    else:

        return "ERROR: unrecognized plot type"

    sns.boxplot(data_grouped, x=x, y=y, hue=hue, dodge=dodge,
                showcaps=True, boxprops={'facecolor':'None'}, showfliers=False, linewidth=2, whiskerprops={'linewidth':1},
                ax=dots, legend=False)

    sns.move_legend(dots, loc='upper left', bbox_to_anchor=(1, 1))
    plt.setp(dots.get_legend().get_texts(), fontsize='15') # for legend text
    plt.setp(dots.get_legend().get_title(), fontsize='20')

    plt.xlabel(xlab, fontweight='bold', fontsize='20')
    plt.ylabel(ylab, fontweight='bold', fontsize='20')
    plt.xticks(fontweight='bold', fontsize='15')
    plt.yticks(fontweight='bold', fontsize='15')

    plt.tight_layout()

    plt.savefig('abundance_v1.png', dpi=300)
    plt.close()

### ------------------MAIN------------------ ###

import seaborn as sns

from matplotlib import pyplot as plt
