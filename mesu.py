#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pymnet as pn
import itertools

def mesu(M,s,output_function=print):
    t = dict()
    for ii,nl in enumerate(M.iter_node_layers()):
        t[nl] = ii
    for nl in M.iter_node_layers():
        S = [{elem_layer} for elem_layer in nl]
        VM_extension = set()
        for neighbor in M[nl]:
            if all(t[addition] > t[nl] for addition in _additions(S,neighbor,t)):
                VM_extension.add(neighbor)
        _multilayer_extend_subgraph(M,s,S,VM_extension,t,nl,output_function)

def _multilayer_extend_subgraph(M,s,S,VM_extension,t,nl,output_function):
    if all(len(S[ii])==s[ii] for ii in range(len(s))):
        if pn.nx.is_connected(pn.transforms.get_underlying_graph(pn.subnet(M,*S))):
            output_function(S)
            return
        else:
            return
    elif not all(len(S[ii])<=s[ii] for ii in range(len(s))):
        return
    N = {neighbor for neighbor in _get_S_neighbors(M,S,t) if t[neighbor] > t[nl]}
    # Exclusive neighborhood to not include nodelayers _within_ the span?
    while VM_extension:
        sigma = VM_extension.pop()
        VM_extension_prime = set(VM_extension)
        new_S = [S[ii].union({sigma[ii]}) for ii in range(len(sigma))]
        new_nodelayers = {new_nl for new_nl in itertools.product(*new_S) if new_nl in t}
        for addition in _additions(S,sigma,t):
            for neighbor in M[addition]:
                if neighbor not in new_nodelayers:
                    if all(t[neighbor_addition] > t[nl] and neighbor_addition not in N for neighbor_addition in _additions(new_S,neighbor,t)):
                        VM_extension_prime.add(neighbor)
        _multilayer_extend_subgraph(M,s,new_S,VM_extension_prime,t,nl,output_function)
    return

def _additions(base_S,added_nl,t):
    # return generator with all nodelayers added as result of adding added_nl to base_S
    total_S = [base_S[ii].union({added_nl[ii]}) for ii in range(len(added_nl))]
    init_subgraph_nls = {nl_init for nl_init in itertools.product(*base_S) if nl_init in t}
    return (nl for nl in itertools.product(*total_S) if (nl in t and not nl in init_subgraph_nls))

def _get_S_neighbors(M,S,t):
    subgraph_nls = {nl for nl in itertools.product(*S) if nl in t}
    for nl in subgraph_nls:
        for neighbor in M[nl]:
            if neighbor not in subgraph_nls:
                yield neighbor
