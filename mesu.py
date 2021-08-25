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
        extension = [set() for elem_layer in nl]
        for neighbor in M[nl]:
            if all(t[addition] > t[nl] for addition in _additions(S,neighbor,t)):
                _add_to_extension(extension,neighbor,S)
        _multilayer_extend_subgraph(M,s,S,extension,t,nl,output_function)

def _multilayer_extend_subgraph(M,s,S,extension,t,nl,output_function):
    if all(len(S[ii])==s[ii] for ii in range(len(s))):
        if _valid(M,S):
            output_function(S)
            return
        else:
            return
#    elif not all(len(S[ii])<=s[ii] for ii in range(len(s))):
#        return
    N = [set() for _ in nl]
    for neighbor in _get_S_neighbors(M,S,t):
        if t[neighbor] > t[nl]:
            for ii,l in enumerate(neighbor):
                N[ii].add(l)
    possible_indices = _candidate_extension_indices(extension,S,s)
    while possible_indices:
        chosen_index = possible_indices[0]
        l = extension[chosen_index].pop()
        possible_indices = _candidate_extension_indices(extension,S,s)
        extension_prime = list(extension)
        # TODO: how to add to just one index in S and keep the rest, without altering original?
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

def _add_to_extension(extension,delta,S,N=None):
    if N is None:
        N = [set() for _ in delta]
    for ii,elem_layer in enumerate(delta):
        if elem_layer not in S[ii] and elem_layer not in S[ii]:
            extension[ii].add(elem_layer)

def _get_S_neighbors(M,S,t):
    subgraph_nls = {nl for nl in itertools.product(*S) if nl in t}
    for nl in subgraph_nls:
        for neighbor in M[nl]:
            if neighbor not in subgraph_nls:
                yield neighbor

def _candidate_extension_indices(extension,S,s):
    return [ii for ii,Ei in enumerate(extension) if Ei and len(S[ii])<s[ii]]

def _valid(M,S):
    subnet = pn.subnet(M,*S)
    valid_elem_layers = [set() for _ in S]
    if pn.nx.is_connected(pn.transforms.get_underlying_graph(subnet)):
        for nl in subnet.iter_node_layers():
            for ii,l in enumerate(nl):
                valid_elem_layers[ii].add(l)
        if valid_elem_layers == S:
            return True
        else:
            return False
    else:
        return False
