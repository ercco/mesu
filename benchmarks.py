#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import time
import matplotlib.pyplot as plt
import helpers

@helpers.persistent
def compare_running_times(subnet_sizes = [(2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)],er_params=([10,5,5],0.05)):
        M = helpers.er_multilayer_any_aspects(*er_params)
        result_times = dict()
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

def run_benchmark(subnet_sizes=((2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)),er_params=((30,10,10),0.0005)):
    persistent_file_name = str(subnet_sizes).replace(' ','')+str(er_params).replace(' ','')
    result_times = compare_running_times(subnet_sizes=subnet_sizes,er_params=er_params,persistent_file=persistent_file_name+'.pickle')
    plot_running_times(result_times,persistent_file_name+'.pdf')






















