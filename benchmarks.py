#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import random
import time
import matplotlib.pyplot as plt
import os,pickle

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

def esu_for_multilayer(M,s,output_function=print):
    t = dict()
    for ii,nl in enumerate(M.iter_node_layers()):
        t[nl] = ii
    for nl in M.iter_node_layers():
        S = [{elem_layer} for elem_layer in nl]
        V_subgraph = {nl}
        extension = set()
        for neighbor in M[nl]:
            if t[neighbor] > t[nl]:
                extension.add(neighbor)
        _esu_extend(M,s,S,V_subgraph,extension,t,nl,output_function)

def _esu_extend(M,s,S,V_subgraph,extension,t,nl,output_function,depth=1):
    if all(len(S[ii])==s[ii] for ii in range(len(s))):
        if mesu.pn.nx.is_connected(mesu.pn.transforms.get_underlying_graph(mesu.pn.subnet(M,*S))):
            output_function(S)
        return
    elif any(len(S[ii])>s[ii] for ii in range(len(s))):
        return
    max_nls = 1
    for ii in s:
        max_nls = max_nls*ii
    if len(V_subgraph) > max_nls:
        return
    N = set()
    for neighbor in mesu._get_S_neighbors(M,S,t):
        N.add(neighbor)
    while extension:
        new_nl = extension.pop()
        S_prime = [S[ii].union({new_nl[ii]}) for ii in range(len(nl))]
        new_exclusive_neighbors = set()
        for neighbor in M[new_nl]:
            if neighbor not in N and neighbor not in V_subgraph and t[neighbor] > t[nl]:
                new_exclusive_neighbors.add(neighbor)
        _esu_extend(M,s,S_prime,V_subgraph.union({new_nl}),extension.union(new_exclusive_neighbors),t,nl,output_function,depth=depth+1)

def er_multilayer_any_aspects(l=[10,4,4],p=0.05):
    M = mesu.pn.MultilayerNetwork(aspects=len(l)-1,directed=False,fullyInterconnected=True)
    for ii,layers_in_aspect in enumerate(l):
        for jj in range(layers_in_aspect):
            M.add_layer(layer=jj,aspect=ii)
    nls = list(M.iter_node_layers())
    for ii in range(len(nls)):
        for jj in range(ii+1,len(nls)):
            if random.random() < p:
                M[nls[ii]][nls[jj]] = 1
    return M

@persistent
def compare_running_times(subnet_sizes = [(2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)],er_params=([10,5,5],0.05)):
        M = er_multilayer_any_aspects(*er_params)
        result_times = dict()
        for subnet_size in subnet_sizes:
            resultset_mesu = set()
            resultset_esu = set()
            mesu_start = time.time()
            mesu.mesu(M,subnet_size,lambda S:resultset_mesu.add(tuple(frozenset(e) for e in S)))
            mesu_end = time.time()
            esu_start = time.time()
            esu_for_multilayer(M,subnet_size,lambda S:resultset_esu.add(tuple(frozenset(e) for e in S)))
            esu_end = time.time()
            assert resultset_mesu == resultset_esu
            result_times[subnet_size] = (mesu_end-mesu_start,esu_end-esu_start)
        return result_times

def plot_running_times(running_times,savename):
    # running times as dict (s1,s2,s3,...): (x,y); x = mesu time, y = esu time
    sorted_subnet_sizes = sorted(running_times.keys())
    mesu_series = []
    esu_series = []
    for s in sorted_subnet_sizes:
        mesu_series.append(running_times[s][0])
        esu_series.append(running_times[s][1])
    plt.plot(mesu_series)
    plt.plot(esu_series)
    plt.legend(['MESU','adapted ESU'])
    plt.xticks(ticks=list(range(len(sorted_subnet_sizes))),labels=[str(s) for s in sorted_subnet_sizes],rotation=45)
    plt.ylabel('Running time')
    plt.xlabel('Subgraph size')
    plt.tight_layout(pad=0.1)
    plt.savefig(savename)
    plt.close('all')

def run_benchmark():
    subnet_sizes = ((2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4))
    er_params = ((10,5,5),0.01)
    persistent_file_name = str(subnet_sizes).replace(' ','')+str(er_params).replace(' ','')
    result_times = compare_running_times(subnet_sizes=subnet_sizes,er_params=er_params,persistent_file=persistent_file_name+'.pickle')
    plot_running_times(result_times,persistent_file_name+'.pdf')






















