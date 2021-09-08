#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import os
import pickle
import itertools
import pymnet as pn
import matplotlib.pyplot as plt
import numpy as np

class persistent(object):
    def __init__(self,function):
        self.function=function
    def __call__(self,*args,**kwargs):
        if 'persistent_file' in kwargs:
            persistent_file=kwargs.pop('persistent_file')
            if os.path.exists(persistent_file): #read from the file
                with open(persistent_file,'rb') as f:
                    result = pickle.load(f)
            else: #create it, write it
                result = self.function(*args,**kwargs)
                with open(persistent_file,'wb') as f:
                    pickle.dump(result,f)
            return result
        else:
            return self.function(*args,**kwargs)

def er_multilayer_any_aspects(l=[10,4,4],p=0.05):
    # generate all nodelayers and then add random edges between them
    M = pn.MultilayerNetwork(aspects=len(l)-1,directed=False,fullyInterconnected=True)
    for ii,layers_in_aspect in enumerate(l):
        for jj in range(layers_in_aspect):
            M.add_layer(layer=jj,aspect=ii)
    nls = list(M.iter_node_layers())
    for ii in range(len(nls)):
        for jj in range(ii+1,len(nls)):
            if random.random() < p:
                M[nls[ii]][nls[jj]] = 1
    return M

def er_multilayer_any_aspects_deg_1_or_greater(l=[5,5,5],p=0.05):
    # only add nodelayers which have at least one edge, net not fully interconnected
    M = pn.MultilayerNetwork(aspects=len(l)-1,directed=False,fullyInterconnected=False)
    possible_nls = list(itertools.product(*(list(range(elem_number)) for elem_number in l)))
    for ii in range(len(possible_nls)):
        for jj in range(ii+1,len(possible_nls)):
            if random.random() < p:
                M[possible_nls[ii]][possible_nls[jj]] = 1
    return M

def heatmap(data, row_labels, col_labels, xlabel, ylabel, title='', ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    cbar.outline.set_visible(False)

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_xticklabels(col_labels)
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    #ax.tick_params(top=True, bottom=False,
    #               labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    #ax.spines[:].set_visible(False)
    ax.set_frame_on(False)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
