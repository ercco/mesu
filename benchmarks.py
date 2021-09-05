#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import time
import matplotlib.pyplot as plt
import helpers

@helpers.persistent
def compare_running_times(subnet_sizes = [(2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)],er_params=([5,5,5],0.1),total_p=None):
    M = helpers.er_multilayer_any_aspects_deg_1_or_greater(*er_params)
    results = dict()
    for subnet_size in subnet_sizes:
        if total_p is None:
            p = None
        else:
            p_depth = sum(s_i-1 for s_i in subnet_size) + 1
            p = [total_p**(1.0/p_depth)] * p_depth
        resultset_mesu = set()
        resultset_esu = set()
        mesu_start = time.time()
        mesu.mesu(M,subnet_size,lambda S:resultset_mesu.add(tuple(frozenset(e) for e in S)),p=p)
        mesu_end = time.time()
        esu_start = time.time()
        mesu.augmented_esu(M,subnet_size,lambda S:resultset_esu.add(tuple(frozenset(e) for e in S)),p=p)
        esu_end = time.time()
        if total_p is not None:
            assert resultset_mesu == resultset_esu
        results[subnet_size] = (mesu_end-mesu_start,esu_end-esu_start,len(resultset_mesu),len(resultset_esu))
    return results

def plot_running_times(running_times,savename):
    # running times as dict (s1,s2,s3,...): (x,y); x = mesu time, y = esu time
    sorted_subnet_sizes = sorted(running_times.keys())
    mesu_series = []
    esu_series = []
    for s in sorted_subnet_sizes:
        mesu_series.append(running_times[s][0])
        esu_series.append(running_times[s][1])
    plt.plot(mesu_series,label='MESU',color='#1f77b4')
    plt.plot(esu_series,label='adapted ESU',color='#ff7f0e')
    plt.legend()
    plt.xticks(ticks=list(range(len(sorted_subnet_sizes))),labels=[str(s) for s in sorted_subnet_sizes],rotation=45)
    y_range = plt.gca().get_ylim()
    for ii,subnet_size in enumerate(sorted_subnet_sizes):
        plt.text(ii,0.025*(y_range[1]-y_range[0]),str(running_times[subnet_size][2]),horizontalalignment='center',color='#1f77b4',fontsize='x-small')
        plt.text(ii,0.0*(y_range[1]-y_range[0]),str(running_times[subnet_size][3]),horizontalalignment='center',color='#ff7f0e',fontsize='x-small')
    plt.ylabel('Running time')
    plt.xlabel('Subgraph size')
    plt.tight_layout(pad=0.1)
    plt.savefig(savename)
    plt.close('all')

def run_benchmark(subnet_sizes=((2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)),er_params=((5,5,5),0.1),total_p=None):
    persistent_file_name = str(subnet_sizes).replace(' ','')+'_'+str(er_params).replace(' ','')+'_'+str(total_p)
    result_times = compare_running_times(subnet_sizes=subnet_sizes,er_params=er_params,total_p=total_p,persistent_file=persistent_file_name+'.pickle')
    plot_running_times(result_times,persistent_file_name+'.pdf')

def run_density_sweep():
    for p in [0.001,0.005,0.01,0.05]:
        run_benchmark(er_params=((5,5,5),p))




















