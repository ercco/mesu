#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import unittest
import time
import scipy,scipy.stats
import pymnet
from pymnet import net,models
from pymnet.sampling import reqs,dumb,esu
import random
import mesu
import pickle
import itertools

### model creating functions

def er_multilayer_partially_interconnected(nodes_by_layer,p,seed=None):
    """Create a one-aspect E-R multilayer network with given nodesets for each
    layer and edge probability p.
    
    Parameters
    ----------
    nodes_by_layer : sequence/iterator of sequences/iterators
        A sequence where each element is a sequence of nodes on a layer.
    p : float 0 <= p <= 1
        The probability that an edge exists between a node-layer pair.
    seed : int, str, bytes or bytearray
        Seed for network generation.
        
    Returns
    -------
    The generated network.
    """
    if seed == None:
        random.seed()
    else:
        random.seed(seed)
    network = pymnet.MultilayerNetwork(aspects=1,fullyInterconnected=False)
    for layer,nodelist in enumerate(nodes_by_layer):
        network.add_layer(layer)
        for node in nodelist:
            network.add_node(node=node,layer=layer)
    numberings = dict()
    for index,nodelayer in enumerate(network.iter_node_layers()):
        numberings[nodelayer] = index
    for nodelayer1 in numberings:
        for nodelayer2 in numberings:
            if numberings[nodelayer1] > numberings[nodelayer2] and random.random() < p:
                network[nodelayer1][nodelayer2] = 1
    return network
    
def random_nodelists(poolsize,nodes_per_layer,layers,seed=None):
    """Draw a random sample of nodes without replacement for each layer
    from a pool of specified size.
    
    Parameters
    ----------
    poolsize : int
        Size of the pool to draw nodes from.
    nodes_per_layer : int
        How many nodes are on each layer.
    layers : int
        How many layers should nodes be drawn for.
    seed : int, str, bytes or bytearray
        Seed for random drawing.
        
    Yields
    ------
    A list of a sample of nodes_per_layer nodes without replacement, times layers.
    """
    if seed == None:
        random.seed()
    else:
        random.seed(seed)
    for _ in range(layers):
        yield random.sample(range(poolsize),nodes_per_layer)
        
### tests

class TestSampling(unittest.TestCase):

    def test_basic_nets_1_aspect(self):
        for func in self.functions_to_test:
            with self.subTest(f=func):
                net1 = net.MultilayerNetwork(aspects=1,fullyInterconnected=False)
                net2 = net.MultilayerNetwork(aspects=1,fullyInterconnected=False)
                net2[1,'X'][1,'Y'] = 1
                net2[1,'X'][2,'X'] = 1
                net3 = net.MultilayerNetwork(aspects=1,fullyInterconnected=False)
                net3[1,'X'][1,'Y'] = 1
                net3[1,'X'][3,'X'] = 1
                net3[1,'Y'][1,'Z'] = 1
                net3[1,'Y'][2,'Z'] = 1
                net4 = net.MultilayerNetwork(aspects=1,fullyInterconnected=True)
                net4[1,'X'][1,'Y'] = 1
                net4[1,'X'][2,'X'] = 1
                net5 = models.full_multilayer(2,['X','Y','Z'])
                net6 = net.MultilayerNetwork(aspects=1,fullyInterconnected=True)
                net6[1,'X'][2,'X'] = 1
                net6[1,'X'][1,'Y'] = 1
                net6[1,'Y'][1,'Z'] = 1
                net6[1,'Z'][2,'Z'] = 1
                net6[2,'Z'][2,'Y'] = 1
                net7 = net.MultilayerNetwork(aspects=1,fullyInterconnected=True)
                net7[1,'X'][2,'X'] = 1
                net7[1,'X'][1,'Y'] = 1
                net7[1,'Y'][1,'Z'] = 1
                net7[1,'Z'][2,'Z'] = 1
                net8 = net.MultilayerNetwork(aspects=1,fullyInterconnected=False)
                net8[1,'X'][1,'Y'] = 1
                net8[1,'X'][2,'X'] = 1
                net8.add_node(2,layer='Y')
                resultlist = []
                func(net1,[1,1],lambda S: resultlist.append(tuple(list(x) for x in S)))
                self.assertEqual(resultlist,[],'\n net1')
                resultlist = []
                func(net2,[1,1],lambda S: resultlist.append(tuple(list(x) for x in S)))
                resultlist.sort()
                self.assertEqual(resultlist,[([1],['X']),([1],['Y']),([2],['X'])],'\n net2 [1,1]')
                resultlist = []
                func(net2,[2,1],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                self.assertEqual(resultlist,[([1,2],['X'])],'\n net2 [2,1]')
                resultlist = []
                func(net2,[2,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                self.assertEqual(resultlist,[([1,2,],['X','Y'])],'\n net2 [2,2]')
                resultlist = []
                func(net3,[2,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([1,2],['Y','Z']),([1,3],['X','Y'])],'\n net3 [2,2]')
                resultlist = []
                func(net3,[3,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                self.assertEqual(resultlist,[],'\n net3 [3,2]')
                resultlist = []
                func(net3,[1,3],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([1],['X','Y','Z'])],'\n net3 [1,3]')
                resultlist = []
                func(net3,[1,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([1],['X','Y']),([1],['Y','Z'])],'\n net3 [1,2]')
                resultlist = []
                func(net3,[3,3],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([1,2,3],['X','Y','Z'])],'\n net3 [3,3]')
                resultlist = []
                func(net4,[2,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[],'\n net4 [2,2]')
                resultlist = []
                func(net5,[2,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([0,1],['X','Y']),([0,1],['X','Z']),([0,1],['Y','Z'])],'\n net5 [2,2]')
                resultlist = []
                func(net5,[2,3],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([0,1],['X','Y','Z'])],'\n net5 [2,3]')
                resultlist = []
                func(net6,[2,3],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([1,2],['X','Y','Z'])],'\n net6 [2,3]')
                resultlist = []
                func(net6,[1,3],lambda S: resultlist.append(tuple(list(x) for x in S)))
                for result in resultlist:
                    result[0].sort()
                    result[1].sort()
                resultlist.sort()
                self.assertEqual(resultlist,[([1],['X','Y','Z'])],'\n net6 [1,3]')
                resultlist = []
                func(net7,[2,3],lambda S: resultlist.append(tuple(list(x) for x in S)))
                self.assertEqual(resultlist,[],'\n net7 [2,3]')
                resultlist = []
                func(net8,[2,2],lambda S: resultlist.append(tuple(list(x) for x in S)))
                self.assertEqual(resultlist,[],'\n net8 [2,2]')

    def test_random_nets_1_aspect(self):
        subnet_sizes = [(2,1),(2,2),(3,2),(3,3)]
        for subnet_size in subnet_sizes:
            for _ in range(10):
                network = er_multilayer_partially_interconnected(random_nodelists(5,5,5),0.4)
                resultlist_dumb = []
                dumb.dumb_enumeration(network,resultlist_dumb,nnodes=subnet_size[0],nlayers=subnet_size[1])
                for result in resultlist_dumb:
                    result[0].sort()
                    result[1].sort()
                resultlist_dumb.sort()
                for func in self.functions_to_test:
                    with self.subTest(f=func):
                        resultlist_esu = []
                        func(network,subnet_size,lambda S: resultlist_esu.append(tuple(list(x) for x in S)))
                        for result in resultlist_esu:
                            result[0].sort()
                            result[1].sort()
                        resultlist_esu.sort()
                        try:
                            self.assertEqual(resultlist_dumb,resultlist_esu)
                        except AssertionError:
                            savelist = [network,resultlist_dumb,resultlist_esu]
                            with open('test_random_nets_1_aspect_fail_'+func.__name__+'.pickle','wb') as f:
                                pickle.dump(savelist,f)
                            raise

    def test_small_nets_exhaustive_1_aspect(self):
        def get_all_subnets(nl_list):
            M = pymnet.MultilayerNetwork(aspects=1,fullyInterconnected=False)
            for ii,nl1 in enumerate(nl_list):
                for nl2 in nl_list[ii+1:]:
                    M[nl1[0],nl1[1]][nl2[0],nl2[1]] = 1
            for m in pymnet.transforms.subnet_iter(M):
                yield m
        def get_all_nl_lists(nnodes,nlayers):
            # powerset with >= 2 elems
            nls = []
            for n in range(nnodes):
                for l in range(nlayers):
                    nls.append((n,l))
            return itertools.chain.from_iterable(itertools.combinations(nls,r) for r in range(2,len(nls)+1))
        def check_all_subnets(nl_list):
            for M in get_all_subnets(nl_list):
                for n in range(1,len(list(M.iter_nodes()))):
                    for l in range(1,len(list(M.iter_layers()))):
                        res_dumb = []
                        dumb.dumb_enumeration(M,res_dumb,nnodes=n,nlayers=l)
                        for result1 in res_dumb:
                            result1[0].sort()
                            result1[1].sort()
                        res_dumb.sort()
                        for func in self.functions_to_test:
                            with self.subTest(f=func):
                                res_esu = []
                                func(M,(n,l),lambda S: res_esu.append(tuple(list(x) for x in S)))
                                for result2 in res_esu:
                                    result2[0].sort()
                                    result2[1].sort()
                                res_esu.sort()
                                try:
                                    self.assertEqual(res_dumb,res_esu)
                                except AssertionError:
                                    savelist = [M,res_dumb,res_esu]
                                    with open('test_small_nets_exhaustive_1_aspect_fail_'+func.__name__+'.pickle','wb') as f:
                                        pickle.dump(savelist,f)
                                    raise
        def run_all_nl_lists(nnodes,nlayers):
            nl_lists = list(get_all_nl_lists(nnodes,nlayers))
            for nl_list in nl_lists:
                check_all_subnets(nl_list)
        run_all_nl_lists(2,2)
                
    def test_esu_insane(self):
        # PyPy recommended for speed
        reqlist = [([1,1],[0]),([1,1],[1]),([1,2],[0]),([1,2],[1]),([1,3],[0]),([1,3],[1]),([2,3],[0]),([2,3],[1]),([2,3],[2]),([3,3],[0]),([3,3],[1]),([3,3],[2]),([3,3],[3])]    
        reqlist = reqlist + [([1,1,1],[0,0,0,0]),([1,1,1],[1,0,0,0]),([1,1,1],[1,1,1,1])]
        reqlist = reqlist + [([2,1,1],[0,0,0,0]),([2,1,1],[1,0,0,0]),([2,1,1],[1,1,1,1])]
        reqlist = reqlist + [([2,2,1],[0,0,0,0]),([2,2,1],[1,0,0,0]),([2,2,1],[2,0,0,0]),([2,2,1],[1,1,0,0]),([2,2,1],[1,0,1,0]),([2,2,1],[1,1,1,1]),([2,2,1],[2,0,0,0]),([2,2,1],[2,1,1,1])]
        for requirement in reqlist:
            for _ in range(100):
                network = er_multilayer_partially_interconnected(random_nodelists(30,10,5),0.05)
                resultlist_dumb = []
                resultlist_esu = []
                dumb.dumb_enumeration(network,resultlist_dumb,sizes=requirement[0],intersections=requirement[1])
                esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=requirement[0],intersections=requirement[1])
                for result in resultlist_dumb:
                    result[0].sort()
                    result[1].sort()
                resultlist_dumb.sort()
                for result in resultlist_esu:
                    result[0].sort()
                    result[1].sort()
                resultlist_esu.sort()
                self.assertEqual(resultlist_dumb,resultlist_esu)
                
    def test_esu_performance(self):
        reqlist = [([1,1],[0]),([1,1],[1]),([1,2],[0]),([1,2],[1]),([1,3],[0]),([1,3],[1]),([2,3],[0]),([2,3],[1]),([2,3],[2]),([3,3],[0]),([3,3],[1]),([3,3],[2]),([3,3],[3])]    
        reqlist = reqlist + [([1,1,1],[0,0,0,0]),([1,1,1],[1,0,0,0]),([1,1,1],[1,1,1,1])]
        reqlist = reqlist + [([2,1,1],[0,0,0,0]),([2,1,1],[1,0,0,0]),([2,1,1],[1,1,1,1])]
        reqlist = reqlist + [([2,2,1],[0,0,0,0]),([2,2,1],[1,0,0,0]),([2,2,1],[2,0,0,0]),([2,2,1],[1,1,0,0]),([2,2,1],[1,0,1,0]),([2,2,1],[1,1,1,1]),([2,2,1],[2,0,0,0]),([2,2,1],[2,1,1,1])]
        network = net.MultilayerNetwork(aspects=1,fullyInterconnected=False)
        nodelayerlist = [(5, 0),(11, 0),(17, 0),(18, 0),(22, 0),
                        (23, 0),(25, 0),(27, 0),(28, 0),(29, 0),
                        (2, 1),(7, 1),(9, 1),(13, 1),(15, 1),
                        (18, 1),(19, 1),(22, 1),(25, 1),(29, 1),
                        (0, 2),(1, 2),(7, 2),(9, 2),(10, 2),
                        (15, 2),(20, 2),(22, 2),(25, 2),(26, 2),
                        (30, 2),(31, 2),(1, 3),(2, 3),(5, 3),
                        (8, 3),(14, 3),(17, 3),(19, 3),(20, 3),
                        (23, 3),(29, 3),(30, 3),(31, 3),(0, 4),
                        (1, 4),(7, 4),(8, 4),(14, 4),(15, 4),
                        (16, 4),(18, 4),(24, 4),(29, 4),(31, 4)]
        edgelist = [(0, 25, 2, 2, 1),(0, 22, 2, 0, 1),(0, 2, 2, 1, 1),(0, 14, 4, 4, 1),(0, 2, 4, 1, 1),(0, 15, 4, 1, 1),(1, 1, 2, 4, 1),
                    (1, 24, 2, 4, 1),(1, 22, 3, 2, 1),(1, 20, 3, 3, 1),(1, 7, 3, 1, 1),(1, 25, 3, 0, 1),(1, 17, 3, 3, 1),(1, 23, 3, 0, 1),
                    (1, 10, 3, 2, 1),(1, 17, 4, 3, 1),(1, 26, 4, 2, 1),(1, 19, 4, 3, 1),(1, 29, 4, 3, 1),(2, 22, 1, 2, 1),(2, 14, 1, 4, 1),
                    (2, 24, 1, 4, 1),(2, 19, 1, 3, 1),(2, 18, 1, 1, 1),(2, 29, 3, 1, 1),(5, 10, 0, 2, 1),(5, 28, 0, 0, 1),(5, 17, 0, 3, 1),
                    (5, 18, 0, 4, 1),(5, 15, 0, 1, 1),(5, 20, 3, 2, 1),(5, 17, 3, 0, 1),(5, 23, 3, 0, 1),(5, 22, 3, 0, 1),(7, 25, 1, 0, 1),
                    (7, 26, 1, 2, 1),(7, 20, 1, 3, 1),(7, 8, 1, 4, 1),(7, 28, 2, 0, 1),(7, 19, 2, 1, 1),(7, 25, 2, 1, 1),(7, 29, 2, 0, 1),
                    (7, 16, 2, 4, 1),(7, 8, 4, 3, 1),(7, 29, 4, 4, 1),(7, 25, 4, 1, 1),(7, 15, 4, 1, 1),(8, 29, 3, 0, 1),(8, 25, 3, 0, 1),
                    (8, 14, 3, 3, 1),(8, 8, 3, 4, 1),(8, 10, 3, 2, 1),(8, 18, 4, 4, 1),(8, 9, 4, 1, 1),(8, 29, 4, 4, 1),(9, 29, 1, 4, 1),
                    (9, 19, 1, 3, 1),(9, 18, 1, 1, 1),(9, 29, 1, 0, 1),(9, 29, 2, 1, 1),(9, 29, 2, 4, 1),(9, 29, 2, 3, 1),(10, 23, 2, 0, 1),
                    (11, 29, 0, 1, 1),(13, 24, 1, 4, 1),(14, 18, 3, 0, 1),(14, 17, 3, 3, 1),(14, 22, 3, 2, 1),(14, 23, 4, 3, 1),(15, 27, 4, 0, 1),
                    (15, 17, 4, 0, 1),(16, 25, 4, 2, 1),(16, 23, 4, 0, 1),(16, 25, 4, 1, 1),(16, 26, 4, 2, 1),(17, 29, 3, 1, 1),(18, 25, 0, 0, 1),
                    (18, 18, 0, 4, 1),(18, 26, 1, 2, 1),(18, 24, 1, 4, 1),(18, 26, 4, 2, 1),(18, 23, 4, 3, 1),(18, 23, 4, 0, 1),(18, 22, 4, 1, 1),
                    (19, 29, 1, 1, 1),(20, 27, 2, 0, 1),(20, 29, 2, 3, 1),(20, 27, 3, 0, 1),(22, 24, 0, 4, 1),(22, 29, 0, 0, 1),(22, 29, 0, 3, 1),
                    (22, 25, 2, 2, 1),(23, 26, 3, 2, 1),(25, 26, 1, 2, 1),(25, 25, 1, 2, 1),(25, 29, 1, 4, 1),(25, 26, 2, 2, 1),(27, 29, 0, 0, 1),
                    (30, 31, 2, 2, 1),(30, 31, 3, 4, 1),(30, 31, 3, 3, 1),(31, 31, 2, 3, 1)]
        for nodelayer in nodelayerlist:
            network.add_node(node=nodelayer[0],layer=nodelayer[1])
        for edge in edgelist:
            network[edge[0],edge[1],edge[2],edge[3]] = 1
        start = time.time()
        for _ in range(10):
            for requirement in reqlist:
                resultlist_esu = []
                esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=requirement[0],intersections=requirement[1])
        print("Time taken "+str(time.time()-start)+" s")
        
    def _statistical_sample(self,network,iterations,motif,p,all_subgraphs):
        samplings = dict()
        for subgraph in all_subgraphs:
            subgraph[0].sort()
            subgraph[1].sort()
            samplings[tuple((tuple(subgraph[0]),tuple(subgraph[1])))] = []
        start = time.time()
        for _ in range(iterations):
            resultlist = []
            esu.sample_multilayer_subgraphs_esu(network,resultlist,p=p,sizes=motif[0],intersections=motif[1])
            for result in resultlist:
                result[0].sort()
                result[1].sort()
            resultlist = [tuple((tuple(result[0]),tuple(result[1]))) for result in resultlist]
            resultlist = set(resultlist)
            for entry in samplings:
                if entry in resultlist:
                    samplings[entry].append(1)
                else:
                    samplings[entry].append(0)
        print("Iterated in: "+str(time.time()-start)+" s")
        return samplings
        
    def test_esu_distribution_width(self,network=None,threshold=0.05,iterations=10000,motif=([2,1],[1]),splitlen=200,p=None,all_subgraphs=None):
        """
        A crude test for checking that the width of the sampling distribution corresponds to the 
        width of the binomial distribution from which the samples should originate. Does a repeated
        sampling of a network and calculates Pr for each instance of a motif in a network:
        
        Pr = the probability that, assuming that the sample is from the corresponding binomial distribution,
        there is at least one value in the sample at least as far away from the expected value as the farthest value
        found in the sample.
        
        This is a crude measure of sampling distribution width relative to the width of the binomial distribution
        that the algorithm should be sampling from if everything is correct. If the algorithm samples a distribution wider
        than the binomial distribution, Pr will be small.
        The check is done for each instance of a motif found in the network specified in the code. The actual default network
        may vary between different compilers/interpreters. Since no multiple correction is used and since crossing the
        threshold doesn't automatically mean that the algorithm isn't working correctly, the test is passed whether
        the threshold is crossed or not. If there are multiple Pr's smaller than a reasonable threshold, this might
        indicate that something is wrong with the algorithm.
        PyPy recommended for speed.
        """
        if network == None:
            network = er_multilayer_partially_interconnected(random_nodelists(100,30,10,seed=1),0.05,seed=1)       
        if p == None:
            req_nodelist_len,req_layerlist_len = reqs.default_calculate_required_lengths(sizes=motif[0],intersections=motif[1])
            p = [0.5] * (req_nodelist_len-1 + req_layerlist_len-1 + 1)
        if all_subgraphs == None:
            all_subgraphs = []
            esu.sample_multilayer_subgraphs_esu(network,all_subgraphs,sizes=motif[0],intersections=motif[1])
        data = self._statistical_sample(network,iterations,motif,p,all_subgraphs)
        outlier_count = 0
        motif_count = len(data)
        for motif_instance in data:
            splitdata = [sum(split) for split in [data[motif_instance][i:i+splitlen] for i in range(0,len(data[motif_instance]),splitlen)]]
            number_of_groups = len(splitdata)
            expected = float(splitlen*scipy.prod(p))
            d_max = max([abs(datapoint-expected) for datapoint in splitdata])
            if 1-abs(scipy.stats.binom.cdf(expected+d_max-1,splitlen,scipy.prod(p)) - scipy.stats.binom.cdf(expected-d_max,splitlen,scipy.prod(p)))**number_of_groups < threshold/float(motif_count):
                outlier_count += 1
        if outlier_count == 0:
            print('No outliers detected at threshold FWER <= '+str(threshold)+', Bonferroni correction used ('+str(motif_count)+' tests, '+str(len(splitdata))+' samples of '+str(splitlen)+' runs each).')
        elif outlier_count == 1:
            print('1 possible outlier at threshold FWER <= '+str(threshold)+', Bonferroni correction used ('+str(motif_count)+' tests, '+str(len(splitdata))+' samples of '+str(splitlen)+' runs each).')
        else:
            print(str(outlier_count)+' possible outliers at threshold FWER <= '+str(threshold)+', Bonferroni correction used ('+str(motif_count)+' tests, '+str(len(splitdata))+' samples of '+str(splitlen)+' runs each).')
        
    def test_different_parameter_sets(self):
        for _ in range(50):
            network = er_multilayer_partially_interconnected(random_nodelists(45,15,5),0.05)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,sizes=[2,1],intersections=[1])
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=[2,1],intersections=[1])
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,sizes=[3,2,2],intersections=1,nnodes=4)
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=[3,2,2],intersections=1,nnodes=4)
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,sizes=[3,2,2],intersections=2,nnodes=4,intersection_type="less_or_equal")
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=[3,2,2],intersections=2,nnodes=4,intersection_type="less_or_equal")
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,sizes=[2,2,2],intersections=[2,2,2,2],nnodes=4,intersection_type="less_or_equal")
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=[2,2,2],intersections=[2,2,2,2],nnodes=4,intersection_type="less_or_equal")
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,nnodes=3,nlayers=3)
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,nnodes=3,nlayers=3)
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,sizes=[2,3,2],intersections=[2,1,None,None],nnodes=4)
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=[2,3,2],intersections=[2,1,None,None],nnodes=4)
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)
            resultlist_dumb = []
            resultlist_esu = []
            dumb.dumb_enumeration(network,resultlist_dumb,sizes=[2,3,2],intersections=[2,1,None,None],nnodes=4,intersection_type="less_or_equal")
            esu.sample_multilayer_subgraphs_esu(network,resultlist_esu,sizes=[2,3,2],intersections=[2,1,None,None],nnodes=4,intersection_type="less_or_equal")
            for result in resultlist_dumb:
                result[0].sort()
                result[1].sort()
            resultlist_dumb.sort()
            for result in resultlist_esu:
                result[0].sort()
                result[1].sort()
            resultlist_esu.sort()
            self.assertEqual(resultlist_dumb,resultlist_esu)

def makesuite(random_nets=True,exhaustive=True):
    suite = unittest.TestSuite()
    TestSampling.functions_to_test = [mesu.mesu,mesu.augmented_esu]
    suite.addTest(TestSampling("test_basic_nets_1_aspect"))
    if random_nets:
        suite.addTest(TestSampling("test_random_nets_1_aspect"))
    if exhaustive:
        suite.addTest(TestSampling("test_small_nets_exhaustive_1_aspect"))
    '''
    if insane:
        suite.addTest(TestSampling("test_esu_insane"))
    if performance:
        suite.addTest(TestSampling("test_esu_performance"))
    if distribution_width:
        suite.addTest(TestSampling("test_esu_distribution_width"))
    if parameter_sets:
        suite.addTest(TestSampling("test_different_parameter_sets"))
    '''
    return suite

def test_sampling(**kwargs):
    suite=makesuite(**kwargs)
    return unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()

if __name__ == '__main__':
    sys.exit(not test_sampling(random_nets=True,exhaustive=True))


