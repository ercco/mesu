#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import os
import pickle
import itertools
import pymnet as pn
import networkx
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import math
import collections

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

def geo_multilayer_any_aspects(l,edges_in_layers,edges_between_layers):
    # TODO: currently does not work with 1 aspect, since layers are ints instead of tuples
    # tempfix: use 1 as number of elementary layers in second aspect
    #
    # l : number of elementary layers in each aspect
    # edges_in_layers : int with approximate number of edges within layers
    # edges_between_layers : int with approximate number of edges between each pair of layers
    M = pn.MultilayerNetwork(aspects=len(l)-1,directed=False,fullyInterconnected=True)
    n = l[0]
    pos = None
    for ii,layers_in_aspect in enumerate(l):
        for jj in range(layers_in_aspect):
            M.add_layer(layer=jj,aspect=ii)
    layers = list(M.iter_layers()) # for getting all pairs later
    for layer in layers:
        netX,pos = get_single_geo_instance(n,edges_in_layers,pos)
        for node in netX.nodes:
            M.add_node(node,layer=layer)
        for e in netX.edges:
            # add edges between nodelayers inside layer
            node1 = e[0]
            node2 = e[1]
            medge = [node1, node2, *layer]
            M[tuple(medge)] = 1
    for layer_pair in itertools.combinations(layers, 2):
        netX,pos = get_single_geo_instance(n,edges_between_layers,pos)
        for e in netX.edges:
            # add edges between nodelayers between layers
            iledge = []
            node1 = e[0]
            node2 = e[1]
            iledge.append(node1)
            iledge.append(node2)
            for ii in range(len(layer_pair[0])):
                iledge.append(layer_pair[0][ii])
                iledge.append(layer_pair[1][ii])
            M[tuple(iledge)] = 1
    return M

def geo_er_mixed_multilayer_any_aspects(l,edges_in_layers,edges_between_layers,p_er_intra,p_er_inter):
    # geo model with er model on top (each layer and each pair of layers)
    #
    # TODO: currently does not work with 1 aspect, since layers are ints instead of tuples
    # tempfix: use 1 as number of elementary layers in second aspect
    #
    # l : number of elementary layers in each aspect
    # edges_in_layers : int with approximate number of edges within layers
    # edges_between_layers : int with approximate number of edges between each pair of layers
    M = pn.MultilayerNetwork(aspects=len(l)-1,directed=False,fullyInterconnected=True)
    n = l[0]
    pos = None
    for ii,layers_in_aspect in enumerate(l):
        for jj in range(layers_in_aspect):
            M.add_layer(layer=jj,aspect=ii)
    layers = list(M.iter_layers()) # for getting all pairs later
    for layer in layers:
        netX,pos = get_single_geo_instance(n,edges_in_layers,pos)
        for node in netX.nodes:
            M.add_node(node,layer=layer)
        for e in netX.edges:
            # add edges between nodelayers inside layer
            node1 = e[0]
            node2 = e[1]
            medge = [node1, node2, *layer]
            M[tuple(medge)] = 1
        # add additional er edges
        all_nodes = list(netX.nodes)
        for node_pair in itertools.combinations(all_nodes,2):
            if random.random() < p_er_intra:
                meredge = [node_pair[0], node_pair[1], *layer]
                M[tuple(meredge)] = 1
    for layer_pair in itertools.combinations(layers, 2):
        netX,pos = get_single_geo_instance(n,edges_between_layers,pos)
        for e in netX.edges:
            # add edges between nodelayers between layers
            iledge = []
            node1 = e[0]
            node2 = e[1]
            iledge.append(node1)
            iledge.append(node2)
            for ii in range(len(layer_pair[0])):
                iledge.append(layer_pair[0][ii])
                iledge.append(layer_pair[1][ii])
            M[tuple(iledge)] = 1
        # add additional er edges
        # consider all possible edge pairs, i.e.
        # bipartite er where number of edges is n^2
        # also including self-edges
        all_nodes = list(netX.nodes)
        for node1 in all_nodes:
            for node2 in all_nodes:
                if random.random() < p_er_inter:
                    ileredge = []
                    ileredge.append(node1)
                    ileredge.append(node2)
                    for ii in range(len(layer_pair[0])):
                        ileredge.append(layer_pair[0][ii])
                        ileredge.append(layer_pair[1][ii])
                    M[tuple(ileredge)] = 1
    return M

def geo_er_mixed_multilayer_any_aspects_mean_deg_parametrization(l,d_geo_intra,d_geo_inter,d_er_intra,d_er_inter):
    # geo model with er model on top (each layer and each pair of layers)
    #
    # TODO: currently does not work with 1 aspect, since layers are ints instead of tuples
    # tempfix: use 1 as number of elementary layers in second aspect
    #
    # l : number of elementary layers in each aspect
    # edges_in_layers : int with approximate number of edges within layers
    # edges_between_layers : int with approximate number of edges between each pair of layers
    #
    # transform to other parametrization:
    # geo
    edges_in_layers = int(np.floor((d_geo_intra*l[0])/2.0))
    edges_between_layers = int(np.floor((l[0]*d_geo_inter)/(np.product(l[1:])-1)))
    # er
    # testing with 0
    p_er_intra = d_er_intra/(float(l[0]-1))
    p_er_inter = d_er_inter/(float(np.product(l[1:])-1)*l[0])
    #
    M = pn.MultilayerNetwork(aspects=len(l)-1,directed=False,fullyInterconnected=True)
    n = l[0]
    pos = None
    for ii,layers_in_aspect in enumerate(l):
        for jj in range(layers_in_aspect):
            M.add_layer(layer=jj,aspect=ii)
    layers = list(M.iter_layers()) # for getting all pairs later
    for layer in layers:
        netX,pos = get_single_geo_instance(n,edges_in_layers,pos)
        for node in netX.nodes:
            M.add_node(node,layer=layer)
        for e in netX.edges:
            # add edges between nodelayers inside layer
            node1 = e[0]
            node2 = e[1]
            medge = [node1, node2, *layer]
            M[tuple(medge)] = 1
        # add additional er edges
        all_nodes = list(netX.nodes)
        for node_pair in itertools.combinations(all_nodes,2):
            if random.random() < p_er_intra:
                meredge = [node_pair[0], node_pair[1], *layer]
                M[tuple(meredge)] = 1
    for layer_pair in itertools.combinations(layers, 2):
        netX,pos = get_single_geo_instance(n,edges_between_layers,pos)
        for e in netX.edges:
            # add edges between nodelayers between layers
            iledge = []
            node1 = e[0]
            node2 = e[1]
            iledge.append(node1)
            iledge.append(node2)
            for ii in range(len(layer_pair[0])):
                iledge.append(layer_pair[0][ii])
                iledge.append(layer_pair[1][ii])
            M[tuple(iledge)] = 1
        # add additional er edges
        # consider all possible edge pairs, i.e.
        # bipartite er where number of edges is n^2
        # also including self-edges
        all_nodes = list(netX.nodes)
        for node1 in all_nodes:
            for node2 in all_nodes:
                if random.random() < p_er_inter:
                    ileredge = []
                    ileredge.append(node1)
                    ileredge.append(node2)
                    for ii in range(len(layer_pair[0])):
                        ileredge.append(layer_pair[0][ii])
                        ileredge.append(layer_pair[1][ii])
                    M[tuple(ileredge)] = 1
    return M

def get_single_geo_instance(n,n_edges,pos):
    # slightly overgenerates the number of edges
    r = math.sqrt(2.2 * float(n_edges) / ((n - 1.0) * n) / math.pi)
    netX = networkx.soft_random_geometric_graph(n, r, pos=pos)
    pos = networkx.get_node_attributes(netX, 'pos')
    return netX, pos

def get_edge_information(M):
    edge_info = dict()
    for layer in M.iter_layers():
        if isinstance(layer,int):
            edge_info[(layer,)] = 0
        else:
            edge_info[layer] = 0
    for layer_pair in itertools.combinations(M.iter_layers(),2):
        unordered_pair = frozenset(layer_pair)
        edge_info[unordered_pair] = 0
    for e in M.edges:
        layer_part = e[2:-1]
        layer1 = tuple(layer_part[0::2])
        layer2 = tuple(layer_part[1::2])
        if layer1 == layer2:
            edge_info[layer1] = edge_info[layer1] + 1
        else:
            edge_info[frozenset((layer1,layer2))] = edge_info[frozenset((layer1,layer2))] + 1
    return edge_info

def load_edgelist(fname):
    M = pn.MultiplexNetwork(couplings='categorical',directed=False,fullyInterconnected=False)
    with open(fname,'r') as f:
        for line in f:
            layer,node1,node2,weight = line.split()
            if node1 != node2:
                M[node1,layer][node2,layer] = 1
    return M

def save_edgelist_cpp_format(M,savename,include_isolated_nls=False):
    # pymnet edge notation: a0, b0, a1, b1, ..., ak, bk
    # where aN is elem layer of node-layer a for aspect N
    # saving notation: a0 a1 a2 ... b0 b1 b2 ...
    # if include_isolated_nls:
    # all edges are written first and lone nodelayers after them (one per line)
    all_nodelayers_with_edges = set()
    with open(savename,'w') as f:
        for e in M.edges:
            nl1 = e[0:(M.aspects+1)*2:2]
            nl2 = e[1:(M.aspects+1)*2:2]
            all_nodelayers_with_edges.add(nl1)
            all_nodelayers_with_edges.add(nl2)
            f.write(' '.join([str(x) for x in nl1+nl2])+'\n')
        if include_isolated_nls:
            for nl in M.iter_node_layers():
                if nl not in all_nodelayers_with_edges:
                    f.write(' '.join([str(y) for y in nl])+'\n')

def load_edgelist_in_cpp_format(savename):
    # only works when first line is an edge (has two nodelayers)
    aspect_flag = False
    with open(savename,'r') as f:
        for line in f:
            line_list = [int(x) for x in line.strip().split()]
            if not aspect_flag:
                n_aspects = int(len(line_list)/2)-1
                M = pn.MultilayerNetwork(aspects=n_aspects,fullyInterconnected=False)
                aspect_flag = True
            nl1 = tuple(line_list[:n_aspects+1])
            nl2 = tuple(line_list[n_aspects+1:])
            if nl2:
                M[nl1][nl2] = 1
            else:
                M.add_node(nl1[0],layer=nl1[1:])
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

def reformat_ppi_file(inputfile,outputfile):
    with open(inputfile,'r') as f:
        with open(outputfile,'w') as g:
            nodes_by_layer = collections.defaultdict(set)
            for line in f:
                layer, node1, node2, edgeweight = line.split()
                nodes_by_layer[layer] = nodes_by_layer[layer].union((node1,node2))
                reformatted_line = ' '.join((node1, layer, node2, layer))
                g.write(reformatted_line+'\n')
            # interlayer edges
            for layer1,layer2 in itertools.combinations(nodes_by_layer.keys(),2):
                for node in nodes_by_layer[layer1].intersection(nodes_by_layer[layer2]):
                    reformatted_line = ' '.join((node,layer1,node,layer2))
                    g.write(reformatted_line+'\n')
