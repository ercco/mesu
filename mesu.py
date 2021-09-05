#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pymnet as pn
import itertools
import copy
import random

##### Nodelayer-based augmented ESU #####

def augmented_esu(M,s,output_function=print,p=None):
    if p is None:
        p = [1] * (sum(s_i-1 for s_i in s) + 1)
    t = dict()
    for ii,nl in enumerate(M.iter_node_layers()):
        t[nl] = ii
    for nl in M.iter_node_layers():
        if random.random() < p[0]:
            S = [{elem_layer} for elem_layer in nl]
            V_subgraph = {nl}
            extension = set()
            for neighbor in M[nl]:
                if t[neighbor] > t[nl]:
                    extension.add(neighbor)
            _augmented_esu_extend(M,s,S,V_subgraph,extension,t,nl,output_function,p)

def _augmented_esu_extend(M,s,S,V_subgraph,extension,t,nl,output_function,p):
    if all(len(S[ii])==s[ii] for ii in range(len(s))):
        if _valid_esu(M,S,V_subgraph,extension,t,nl):
            output_function(S)
        return
    N = set()
    for neighbor in _get_V_subgraph_neighbors(M,V_subgraph):
        N.add(neighbor)
    while extension:
        new_nl = extension.pop()
        S_prime = [S[ii].union({new_nl[ii]}) for ii in range(len(nl))]
        # check if we go over any aspect
        if any(len(S_prime[ii])>s[ii] for ii in range(len(s))):
            continue
        # find probability to proceed, if S is not increased at all then always proceed (prob=1)
        first_p_index = sum(len(S_i)-1 for S_i in S)+1
        prob = 1
        for ii in range(first_p_index,first_p_index+sum(len(S_prime[jj])-len(S[jj]) for jj in range(len(S)))):
            prob = prob*p[ii]
        if random.random() < prob:
            new_exclusive_neighbors = set()
            for neighbor in M[new_nl]:
                if neighbor not in N and neighbor not in V_subgraph and t[neighbor] > t[nl]:
                    new_exclusive_neighbors.add(neighbor)
            _augmented_esu_extend(M,s,S_prime,V_subgraph.union({new_nl}),extension.union(new_exclusive_neighbors),t,nl,output_function,p)

##### Aspect-based MESU #####

def mesu(M,s,output_function=print,p=None):
    if p is None:
        p = [1] * (sum(s_i-1 for s_i in s) + 1)
    t = dict()
    for ii,nl in enumerate(M.iter_node_layers()):
        t[nl] = ii
    for nl in M.iter_node_layers():
        if random.random() < p[0]:
            S = [{elem_layer} for elem_layer in nl]
            extension = [set() for elem_layer in nl]
            for neighbor in M[nl]:
                if all(t[addition] > t[nl] for addition in _additions(S,neighbor,t)):
                    _add_to_extension(extension,neighbor,S)
            _multilayer_extend_subgraph(M,s,S,extension,t,nl,output_function,p)

def _multilayer_extend_subgraph(M,s,S,extension,t,nl,output_function,p):
    if all(len(S[ii])==s[ii] for ii in range(len(s))):
        if _valid(M,S):
            output_function(S)
        return
    N = [set() for _ in nl]
    for neighbor in _get_S_neighbors(M,S,t):
        if all(t[addition] > t[nl] for addition in _additions(S,neighbor,t)):
            for ii,l in enumerate(neighbor):
                N[ii].add(l)
    possible_indices = _candidate_extension_indices(extension,S,s)
    while possible_indices:
        chosen_index = possible_indices[0]
        l = extension[chosen_index].pop()
        possible_indices = _candidate_extension_indices(extension,S,s)
        if random.random() < p[sum(len(S_i)-1 for S_i in S)+1]:
            extension_prime = copy.deepcopy(extension)
            S_prime = copy.deepcopy(S)
            S_prime[chosen_index].add(l)
            new_nodelayers = {new_nl for new_nl in itertools.product(*S_prime) if new_nl in t}
            dummy_nl = tuple(l if ii==chosen_index else None for ii in range(len(nl)))
            for addition in _additions(S,dummy_nl,t):
                for neighbor in M[addition]:
                    if neighbor not in new_nodelayers:
                        if all(t[neighbor_addition] > t[nl] for neighbor_addition in _additions(S_prime,neighbor,t)):
                            _add_to_extension(extension_prime,neighbor,S_prime,N)
            _multilayer_extend_subgraph(M,s,S_prime,extension_prime,t,nl,output_function,p)
    return

##### Helpers #####

def _additions(base_S,added_nl,t):
    # return generator with all nodelayers added as result of adding added_nl to base_S
    # to add only a portion of aspects in added_nl, set the rest to None
    total_S = [base_S[ii].union({added_nl[ii]}) if added_nl[ii] is not None else base_S[ii] for ii in range(len(added_nl))]
    init_subgraph_nls = {nl_init for nl_init in itertools.product(*base_S) if nl_init in t}
    return (nl for nl in itertools.product(*total_S) if (nl in t and not nl in init_subgraph_nls))

def _add_to_extension(extension,delta,S,N=None):
    if N is None:
        N = [set() for _ in delta]
    for ii,elem_layer in enumerate(delta):
        if elem_layer not in S[ii] and elem_layer not in N[ii]:
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

def _get_V_subgraph_neighbors(M,V_subgraph):
    for nl in V_subgraph:
        for neighbor in M[nl]:
            if neighbor not in V_subgraph:
                yield neighbor

def _valid_esu(M,S,V_subgraph,extension,t,nl):
    subnet = pn.subnet(M,*S)
    if pn.nx.is_connected(pn.transforms.get_underlying_graph(subnet)) and all(t[sub_nl] >= t[nl] for sub_nl in subnet.iter_node_layers()):
        for spanning_path_nl in V_subgraph:
            if any(neighbor not in V_subgraph and neighbor not in extension for neighbor in subnet[spanning_path_nl]):
                return False
        return True
    else:
        return False

##### Naive algorithm for comparison #####

def naive(M,s,output_function=print):
    # must check for empty layers since connectivity returns ok for them
    for S in itertools.product(*map(lambda i:list(itertools.combinations(M.slices[i],s[i])),range(M.aspects+1))):
        try:
            subnet = pn.subnet(M,*S)
            if pn.nx.is_connected(pn.transforms.get_underlying_graph(subnet)):
                valid_elem_layers = [set() for _ in S]
                for nl in subnet.iter_node_layers():
                    for ii,l in enumerate(nl):
                        valid_elem_layers[ii].add(l)
                if valid_elem_layers == [set(S[jj]) for jj in range(len(S))]:
                    output_function(S)
        except pn.nx.networkx.NetworkXPointlessConcept:
            pass




































