#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pymnet as pn
import itertools

def mesu(M,s):
    t = dict()
    for ii,nl in enumerate(M.iter_node_layers()):
        t[nl] = ii
    for nl in M.iter_node_layers():
        S = [{elem_layer} for elem_layer in nl]
        VM_extension = set()
        for neighbor in M[nl]:
            if all(t[addition] > t[nl] for addition in _additions(S,neighbor,t)):
                VM_extension.add(neighbor)
        _multilayer_extend_subgraph(M,s,S,VM_extension,t,nl)
        
def _multilayer_extend_subgraph(M,s,S,VM_extension,t,nl):
    if all(len(S[ii])==s[ii] for ii in range(len(s))):
        if pn.nx.is_connected(pn.transforms.get_underlying_graph(pn.subnet(M,*S))):
            print(S)
            return
        else:
            return
    elif not all(len(S[ii])<=s[ii] for ii in range(len(s))):
        return
    N = {} # TODO
    # Exclusive neighborhood to not include nodelayers _within_ the span?
            
def _additions(base_S,added_nl,t):
    # return generator with all nodelayers added as result of adding added_nl to base_S
    total_S = [base_S[ii].union({added_nl[ii]}) for ii in range(len(added_nl))]
    init_subgraph_nls = {nl_init for nl_init in itertools.product(*base_S) if nl_init in t}
    return (nl for nl in itertools.product(*total_S) if (nl in t and not nl in init_subgraph_nls))