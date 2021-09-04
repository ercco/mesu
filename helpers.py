#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import os
import pickle
import pymnet as pn

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
