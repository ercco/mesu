#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import time
import matplotlib.pyplot as plt
import numpy as np
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
        if total_p is None:
            assert resultset_mesu == resultset_esu
        results[subnet_size] = (mesu_end-mesu_start,esu_end-esu_start,len(resultset_mesu),len(resultset_esu))
    return results

def plot_running_times(running_times,savename,y_log=False,legend=True):
    # running times as dict (s1,s2,s3,...): (x,y,z,w); x = mesu time, y = esu time, z = number of mesu subg's, w = number of esu subg's
    # or as iterable of such dicts, in which case mean and stddev are plotted
    if isinstance(running_times,dict):
        sorted_subnet_sizes = sorted(running_times.keys())
        mesu_series = []
        esu_series = []
        mesu_errors = None
        esu_errors = None
        mesu_numbers = []
        esu_numbers = []
        for s in sorted_subnet_sizes:
            mesu_series.append(running_times[s][0])
            esu_series.append(running_times[s][1])
            mesu_numbers.append(running_times[s][2])
            esu_numbers.append(running_times[s][3])
    else:
        sorted_subnet_sizes = sorted(running_times[0].keys())
        mesu_series = []
        esu_series = []
        mesu_errors = []
        esu_errors = []
        mesu_numbers = []
        esu_numbers = []
        mesu_number_errors = []
        esu_number_errors = []
        for s in sorted_subnet_sizes:
            mesu_iters = ([],[])
            esu_iters = ([],[])
            for running_time_dict in running_times:
                mesu_iters[0].append(running_time_dict[s][0])
                mesu_iters[1].append(running_time_dict[s][2])
                esu_iters[0].append(running_time_dict[s][1])
                esu_iters[1].append(running_time_dict[s][3])
            mesu_series.append(np.mean(mesu_iters[0]))
            esu_series.append(np.mean(esu_iters[0]))
            mesu_errors.append(np.std(mesu_iters[0]))
            esu_errors.append(np.std(esu_iters[0]))
            mesu_numbers.append(np.mean(mesu_iters[1]))
            esu_numbers.append(np.mean(esu_iters[1]))
            mesu_number_errors.append(np.std(mesu_iters[1]))
            esu_number_errors.append(np.std(esu_iters[1]))
    # sort by maximum for each subnet size
    max_arg_list = np.argsort([max(m,e) for m,e in zip(mesu_series,esu_series)])
    mesu_series = [mesu_series[ii] for ii in max_arg_list]
    esu_series = [esu_series[ii] for ii in max_arg_list]
    mesu_numbers = [mesu_numbers[ii] for ii in max_arg_list]
    esu_numbers = [esu_numbers[ii] for ii in max_arg_list]
    if not isinstance(running_times,dict):
        mesu_errors = [mesu_errors[ii] for ii in max_arg_list]
        esu_errors = [esu_errors[ii] for ii in max_arg_list]
        mesu_number_errors = [mesu_number_errors[ii] for ii in max_arg_list]
        esu_number_errors = [esu_number_errors[ii] for ii in max_arg_list]
    sorted_subnet_sizes = [sorted_subnet_sizes[ii] for ii in max_arg_list]
    plt.errorbar(x=range(len(mesu_series)),y=mesu_series,yerr=mesu_errors,label='MESU',color='#1f77b4')
    plt.errorbar(x=range(len(esu_series)),y=esu_series,yerr=esu_errors,label='adapted ESU',color='#ff7f0e')
    if legend:
        plt.legend()
    if y_log:
        plt.gca().set_yscale('log')
    plt.xticks(ticks=list(range(len(sorted_subnet_sizes))),labels=[str(s) for s in sorted_subnet_sizes],rotation=45)
    y_range = plt.gca().get_ylim()
    for ii,subnet_size in enumerate(sorted_subnet_sizes):
        if isinstance(running_times,dict):
            plt.text(ii,0.025*(y_range[1]-y_range[0]),str(mesu_numbers[ii]),horizontalalignment='center',color='#1f77b4',fontsize='x-small')
            plt.text(ii,0.0*(y_range[1]-y_range[0]),str(esu_numbers[ii]),horizontalalignment='center',color='#ff7f0e',fontsize='x-small')
        else:
            plt.text(ii,0.025*(y_range[1]-y_range[0]),str(mesu_numbers[ii])+r'$\pm$'+str(np.around(mesu_number_errors[ii],1)),horizontalalignment='center',color='#1f77b4',fontsize='x-small')
            plt.text(ii,0.0*(y_range[1]-y_range[0]),str(esu_numbers[ii])+r'$\pm$'+str(np.around(esu_number_errors[ii],1)),horizontalalignment='center',color='#ff7f0e',fontsize='x-small')
    plt.ylabel('Running time')
    plt.xlabel('Subgraph size')
    plt.tight_layout(pad=0.1)
    plt.savefig(savename)
    plt.close('all')

def run_benchmark(subnet_sizes=((2,2,2),(2,2,3),(2,3,3),(3,3,3),(3,3,4)),er_params=((5,5,5),0.1),total_p=None,iter_label='',return_result=False):
    if iter_label:
        iter_label = '_'+iter_label
    persistent_file_name = str(subnet_sizes).replace(' ','')+'_'+str(er_params).replace(' ','')+'_'+str(total_p)+iter_label
    result_times = compare_running_times(subnet_sizes=subnet_sizes,er_params=er_params,total_p=total_p,persistent_file=persistent_file_name+'.pickle')
    plot_running_times(result_times,persistent_file_name+'.pdf')
    if return_result:
        return result_times

def run_density_sweep(total_p=None):
    subnet_sizes = ((2,1,1),(2,2,1),(2,2,2),(3,1,1),(3,2,1),(3,2,2))
    if total_p is None:
        for p in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            run_benchmark(subnet_sizes=subnet_sizes,er_params=((5,5,5),p))
    else:
        all_result_times = dict()
        for p in [0.1,0.2,0.3,0.4,0.5]:
            all_result_times[p] = []
            for ii in range(10):
                all_result_times[p].append(run_benchmark(subnet_sizes=subnet_sizes,er_params=((10,10,10),p),total_p=total_p,iter_label=str(ii),return_result=True))
            fig_savename = str(subnet_sizes).replace(' ','')+'_'+str(((10,10,10),p)).replace(' ','')+'_'+str(total_p)+'_alliters.pdf'
            plot_running_times(all_result_times[p],fig_savename,y_log=True)

def run_heatmap_sweep(subnet_size,total_p_10):
    net_sizes = [[10]*len(subnet_size),[20]*len(subnet_size),[30]*len(subnet_size),[40]*len(subnet_size),[50]*len(subnet_size)]
    # avg_deg = 100,150,200,250,300
    densities = [[(100+ii*50)/(np.product(net_size)-1) for ii in range(0,5)] for net_size in net_sizes]
    total_ps = [total_p_10/((net_size[0]/net_sizes[0][0])**len(net_size)) for net_size in net_sizes]
    for ii,net_size in enumerate(net_sizes):
        for density in densities[ii]:
            for iter_label in range(10):
                run_benchmark(subnet_sizes=(subnet_size,),er_params=(net_size,density),total_p=total_ps[ii],iter_label=str(iter_label))


















