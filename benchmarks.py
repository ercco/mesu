#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import time
import matplotlib.pyplot as plt
import helpers

@helpers.persistent
def compare_running_times(subnet_sizes = [(2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)],er_params=([5,5,5],0.1)):
        M = helpers.er_multilayer_any_aspects_deg_1_or_greater(*er_params)
        results = dict()
        for subnet_size in subnet_sizes:
            resultset_mesu = set()
            resultset_esu = set()
            mesu_start = time.time()
            mesu.mesu(M,subnet_size,lambda S:resultset_mesu.add(tuple(frozenset(e) for e in S)))
            mesu_end = time.time()
            esu_start = time.time()
            mesu.augmented_esu(M,subnet_size,lambda S:resultset_esu.add(tuple(frozenset(e) for e in S)))
            esu_end = time.time()
            assert resultset_mesu == resultset_esu
            results[subnet_size] = (mesu_end-mesu_start,esu_end-esu_start,len(resultset_mesu))
        return results

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
    y_range = plt.gca().get_ylim()
    for ii,subnet_size in enumerate(sorted_subnet_sizes):
        plt.text(ii,0.03*(y_range[1]-y_range[0]),str(running_times[subnet_size][2]),horizontalalignment='center',fontsize='x-small')
    plt.ylabel('Running time')
    plt.xlabel('Subgraph size')
    plt.tight_layout(pad=0.1)
    plt.savefig(savename)
    plt.close('all')

def run_benchmark(subnet_sizes=((2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)),er_params=((5,5,5),0.1)):
    persistent_file_name = str(subnet_sizes).replace(' ','')+str(er_params).replace(' ','')
    result_times = compare_running_times(subnet_sizes=subnet_sizes,er_params=er_params,persistent_file=persistent_file_name+'.pickle')
    plot_running_times(result_times,persistent_file_name+'.pdf')

def run_density_sweep():
    for p in [0.1,0.3,0.5,0.7,0.9]:
        run_benchmark(er_params=((5,5,5),p))




















