#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np
import collections
import helpers
import os
from cpp_sanity_check import read_temp_file_res
import pymnet as pn

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

def plot_running_times(running_times,savename,y_log=False,legend=True,plot_to_ax=None,linestyle='solid',show_numbers=True):
    # running times as dict (s1,s2,s3,...): (x,y,z,w); x = mesu time, y = esu time, z = number of mesu subg's, w = number of esu subg's
    # or as iterable of such dicts, in which case mean and stddev are plotted
    if plot_to_ax is None:
        fig,ax = plt.subplots()
    else:
        ax = plot_to_ax
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
    ax.errorbar(x=range(len(mesu_series)),y=mesu_series,yerr=mesu_errors,label='MESU',color='#1f77b4',linestyle=linestyle)
    ax.errorbar(x=range(len(esu_series)),y=esu_series,yerr=esu_errors,label='adapted ESU',color='#ff7f0e',linestyle=linestyle)
    if legend:
        ax.legend()
    if y_log:
        ax.set_yscale('log')
    #plt.xticks(ticks=list(range(len(sorted_subnet_sizes))),labels=[str(s) for s in sorted_subnet_sizes],rotation=45)
    ax.set_xticks(ticks=list(range(len(sorted_subnet_sizes))))
    ax.set_xticklabels([str(s) for s in sorted_subnet_sizes],rotation=45)
    y_range = ax.get_ylim()
    if show_numbers:
        for ii,subnet_size in enumerate(sorted_subnet_sizes):
            if isinstance(running_times,dict):
                ax.text(ii,mesu_series[ii],str(mesu_numbers[ii]),horizontalalignment='center',color='#1f77b4',fontsize='xx-small')
                ax.text(ii,esu_series[ii],str(esu_numbers[ii]),horizontalalignment='center',color='#ff7f0e',fontsize='xx-small')
            else:
                ax.text(ii,mesu_series[ii],str(mesu_numbers[ii])+r'$\pm$'+str(np.around(mesu_number_errors[ii],1)),horizontalalignment='center',color='#1f77b4',fontsize='xx-small')
                ax.text(ii,esu_series[ii],str(esu_numbers[ii])+r'$\pm$'+str(np.around(esu_number_errors[ii],1)),horizontalalignment='center',color='#ff7f0e',fontsize='xx-small')
    if plot_to_ax is None:
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

def run_density_sweep(subnet_sizes=((2,1,1),(2,2,1),(2,2,2),(3,1,1),(3,2,1),(3,2,2)),p_er=[0.1,0.2,0.3,0.4,0.5],total_p=None):
    if total_p is None:
        for p in p_er:
            run_benchmark(subnet_sizes=subnet_sizes,er_params=((5,5,5),p))
    else:
        all_result_times = dict()
        combined_plot_fig,combined_plot_ax = plt.subplots()
        for p in p_er:
            all_result_times[p] = []
            for ii in range(10):
                all_result_times[p].append(run_benchmark(subnet_sizes=subnet_sizes,er_params=((10,10,10),p),total_p=total_p,iter_label=str(ii),return_result=True))
            fig_savename = str(subnet_sizes).replace(' ','')+'_'+str(((10,10,10),p)).replace(' ','')+'_'+str(total_p)+'_alliters.pdf'
            plot_running_times(all_result_times[p],fig_savename,y_log=True)
            if p in (0.1,0.5):
                if p == 0.1:
                    linestyle = 'dashed'
                if p == 0.5:
                    linestyle = 'solid'
                plot_running_times(all_result_times[p],savename='should_not_exist',y_log=True,legend=False,plot_to_ax=combined_plot_ax,linestyle=linestyle)
        combined_plot_ax.set_ylabel('Running time')
        combined_plot_ax.set_xlabel('Subgraph size')
        combined_plot_fig.tight_layout(pad=0.1)
        custom_legend_lines = [Line2D([0],[0],color='#1f77b4'), Line2D([0],[0],color='#ff7f0e'),Line2D([0],[0],color='black',linestyle='dashed'),Line2D([0],[0],color='black',linestyle='solid')]
        combined_plot_ax.legend(custom_legend_lines,['MESU','AESU','density = 0.1','density = 0.5'])
        combined_savename = str(subnet_sizes).replace(' ','')+'_'+str(((10,10,10),(0.1,0.3,0.5))).replace(' ','')+'_'+str(total_p)+'_combined.pdf'
        combined_plot_fig.savefig(combined_savename,bbox_inches='tight')
        plt.close('all')

def plot_combined_with_different_sampling_ps():
    subnet_sizes = ((2,1,1),(2,2,1),(2,2,2),(3,1,1),(3,2,1),(3,2,2),(3,3,2))
    combined_plot_fig,combined_plot_ax = plt.subplots()
    for er_params,total_p in [(((10,10,10),0.1),0.01),(((10,10,10),0.5),0.0001)]:
        alliters = []
        for ii in range(10):
            alliters.append(run_benchmark(subnet_sizes=subnet_sizes,er_params=er_params,total_p=total_p,iter_label=str(ii),return_result=True))
        if er_params[1] == 0.1:
            linestyle = 'dashed'
        if er_params[1] == 0.5:
            linestyle = 'solid'
        plot_running_times(alliters,savename='should_not_exist',y_log=True,legend=False,plot_to_ax=combined_plot_ax,linestyle=linestyle,show_numbers=False)
    combined_plot_ax.set_ylabel('Running time')
    combined_plot_ax.set_xlabel('Subnetwork size')
    combined_plot_fig.tight_layout(pad=0.1)
    custom_legend_lines = [Line2D([0],[0],color='#1f77b4'), Line2D([0],[0],color='#ff7f0e'),Line2D([0],[0],color='black',linestyle='dashed'),Line2D([0],[0],color='black',linestyle='solid')]
    combined_plot_ax.legend(custom_legend_lines,['A-MESU','NL-MESU',r'density = 0.1, $p_{tot}$ = 0.01',r'density = 0.5, $p_{tot}$ = 0.0001'])
    combined_savename = str(subnet_sizes).replace(' ','')+'_'+str(((10,10,10),(0.1,0.5))).replace(' ','')+'_'+str((0.01,0.0001))+'_nonumbers_combined.pdf'
    combined_plot_fig.savefig(combined_savename,bbox_inches='tight')
    plt.close('all')

def run_sampling_prob_sweep(subnet_sizes=((2,1,1),(2,2,1),(2,2,2),(3,1,1),(3,2,1),(3,2,2)),p_er=0.1,total_p=[0.02,0.04,0.06,0.08,0.10]):
    all_result_times = dict()
    for sampling_p in total_p:
        all_result_times[sampling_p] = []
        for ii in range(10):
            all_result_times[sampling_p].append(run_benchmark(subnet_sizes=subnet_sizes,er_params=((10,10,10),p_er),total_p=sampling_p,iter_label=str(ii),return_result=True))
        fig_savename = str(subnet_sizes).replace(' ','')+'_'+str(((10,10,10),p_er)).replace(' ','')+'_'+str(sampling_p)+'_alliters.pdf'
        plot_running_times(all_result_times[sampling_p],fig_savename,y_log=True)

def run_heatmap_sweep(subnet_size,total_p_10):
    net_sizes = [[10]*len(subnet_size),[20]*len(subnet_size),[30]*len(subnet_size),[40]*len(subnet_size),[50]*len(subnet_size)]
    # avg_deg = 100,150,200,250,300
    # densities = [[(100+ii*50)/(np.product(net_size)-1) for ii in range(0,5)] for net_size in net_sizes]
    # INSTEAD: just create nets by simple density?
    densities = [[0.10,0.11,0.12,0.13,0.14,0.15]]*len(net_sizes)
    total_ps = [total_p_10/((net_size[0]/net_sizes[0][0])**len(net_size)) for net_size in net_sizes]
    for ii,net_size in enumerate(net_sizes):
        for density in densities[ii]:
            for iter_label in range(10):
                run_benchmark(subnet_sizes=(subnet_size,),er_params=(net_size,density),total_p=total_ps[ii],iter_label=str(iter_label))

def run_heatmap_sweep_density_and_sampling_prob(subnet_size):
    densities = [0.1,0.15,0.2,0.25,0.3]
    sampling_probs = [0.001,0.004,0.007,0.01]
    resultgrid = []
    for d in densities:
        x_direction_list = []
        for sampling_p in sampling_probs:
            iterlist = []
            for iter_label in range(10):
                iterlist.append(run_benchmark(subnet_sizes=(subnet_size,),er_params=((10,10,10),d),total_p=sampling_p,iter_label=str(iter_label),return_result=True))
            x_direction_list.append(iterlist)
        resultgrid.append(x_direction_list)
    savename = 'heatmap.pdf'
    plot_heatmap_from_iters(resultgrid,sampling_probs,densities,subnet_size,xlabel='Sampling p',ylabel='Density',title='',savename=savename,center=1)
    #savename = str(subnet_size).replace(' ','')+'_'+str(((10,10,10),densities)).replace(' ','')+'_'+str(sampling_probs).replace(' ','')+'_heatmap.pdf'
    #plot_heatmap_from_iters(resultgrid,sampling_probs,densities,subnet_size,xlabel='Sampling p',ylabel='Density',title=str(subnet_size)+', net: '+str((10,10,10)),savename=savename,center=1)

def run_heatmap_sweep_density_and_net_size(subnet_size):
    densities = [0.1,0.15,0.2,0.25,0.3]
    net_sizes = [(10,10,10),(11,11,11),(12,12,12),(13,13,13),(14,14,14)]
    sampling_p = 0.001
    resultgrid = []
    for d in densities:
        x_direction_list = []
        for net_size in net_sizes:
            iterlist = []
            for iter_label in range(10):
                iterlist.append(run_benchmark(subnet_sizes=(subnet_size,),er_params=(net_size,d),total_p=sampling_p,iter_label=str(iter_label),return_result=True))
            x_direction_list.append(iterlist)
        resultgrid.append(x_direction_list)
    savename = str(subnet_size).replace(' ','')+'_'+str(net_sizes).replace(' ','')+'_'+str(sampling_p).replace(' ','')+'_heatmap.pdf'
    plot_heatmap_from_iters(resultgrid,net_sizes,densities,subnet_size,xlabel='Net size',ylabel='Density',title=str(subnet_size)+', sampling p: '+str(sampling_p),savename=savename,center=1)


def plot_heatmap_from_iters(resultgrid,x_axis,y_axis,subnet_size,xlabel,ylabel,title,savename='',center=1):
    heatvaluegrid = []
    for jj,_ in enumerate(resultgrid):
        heat_x_direction_list = []
        for ii,_ in enumerate(resultgrid[jj]):
            mesu_iters = ([],[])
            esu_iters = ([],[])
            for kk,_ in enumerate(resultgrid[jj][ii]):
                mesu_iters[0].append(resultgrid[jj][ii][kk][subnet_size][0])
                mesu_iters[1].append(resultgrid[jj][ii][kk][subnet_size][2])
                esu_iters[0].append(resultgrid[jj][ii][kk][subnet_size][1])
                esu_iters[1].append(resultgrid[jj][ii][kk][subnet_size][3])
            heat_x_direction_list.append(np.mean(mesu_iters[0])/np.mean(esu_iters[0]))
        heatvaluegrid.append(heat_x_direction_list)
    fig,ax = plt.subplots()
    helpers.heatmap(np.array(heatvaluegrid),[str(y) for y in y_axis],[str(x) for x in x_axis],xlabel=xlabel,ylabel=ylabel,title=title,ax=ax,cbarlabel='T(A-MESU)/T(NL-MESU)',norm=colors.TwoSlopeNorm(vcenter=center),cmap='PuOr')
    fig.tight_layout(pad=0.1)
    if savename:
        fig.savefig(savename,bbox_inches='tight')
    else:
        plt.show()

##### Protein-protein interaction nets

@helpers.persistent
def compare_running_times_data(fname,subnet_sizes=[(3,3)],total_p=0.0001):
    M = helpers.load_edgelist(fname)
    net_statistics = (len(list(M.iter_node_layers())),len(list(M.edges)),len(list(M.iter_nodes())),len(list(M.iter_layers())))
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
    return (results,net_statistics)

def run_times_for_example_data():
    subnet_sizes = [(4,3)]
    base_p = 0.00001
    result_times = []
    net_statistics = []
    savename = 'data_'+str(subnet_sizes).replace(' ','')+'_scatter.pdf'
    data = [('multiplex_pp_data/Arabidopsis_Multiplex_Genetic/Dataset/arabidopsis_genetic_multiplex.edges',base_p,'arabidopsis'),
            ('multiplex_pp_data/Bos_Multiplex_Genetic/Dataset/bos_genetic_multiplex.edges',base_p,'bos'),
            ('multiplex_pp_data/Candida_Multiplex_Genetic/Dataset/candida_genetic_multiplex.edges',base_p,'candida'),
            ('multiplex_pp_data/Celegans_Multiplex_Genetic/Dataset/celegans_genetic_multiplex.edges',base_p,'celegans'),
            ('multiplex_pp_data/Drosophila_Multiplex_Genetic/Dataset/drosophila_genetic_multiplex.edges',base_p,'drosophila'),
            ('multiplex_pp_data/Gallus_Multiplex_Genetic/Dataset/gallus_genetic_multiplex.edges',base_p,'gallus'),
            ('multiplex_pp_data/Mus_Multiplex_Genetic/Dataset/mus_genetic_multiplex.edges',base_p,'mus'),
            ('multiplex_pp_data/Plasmodium_Multiplex_Genetic//Dataset/plasmodium_genetic_multiplex.edges',base_p,'plasmodium'),
            ('multiplex_pp_data/Rattus_Multiplex_Genetic/Dataset/rattus_genetic_multiplex.edges',base_p,'rattus'),
            ('multiplex_pp_data/SacchCere_Multiplex_Genetic/Dataset/sacchcere_genetic_multiplex.edges',base_p,'sacchcere'),
            ('multiplex_pp_data/SacchPomb_Multiplex_Genetic/Dataset/sacchpomb_genetic_multiplex.edges',base_p,'sacchpomb')]
    for d in data:
        result_tot = compare_running_times_data(fname=d[0],subnet_sizes=subnet_sizes,total_p=d[1],persistent_file='data_'+d[2]+'_'+str(subnet_sizes[0]).replace(' ','')+'.pickle')
        result_times.append(result_tot[0])
        net_statistics.append(result_tot[1])
        print(d[2])
        print('layers: '+str(result_tot[1][3]))
        print('nodes: '+str(result_tot[1][2]))
        print('nodelayers: '+str(result_tot[1][0]))
    x = [stat[0] for stat in net_statistics]
    y = [stat[1] for stat in net_statistics]
    color = [time[subnet_sizes[0]][0]/time[subnet_sizes[0]][1] for time in result_times]
    sc = plt.scatter(x,y,s=400,c=color,norm=colors.TwoSlopeNorm(vcenter=1.0),cmap='BrBG_r',edgecolors='black')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    cb = plt.colorbar(sc)
    cb.ax.set_ylabel(r'$T_{A-MESU}/T_{NL-MESU}$', rotation=-90, va="bottom")
    plt.xlabel('Number of nodelayers')
    plt.ylabel('Number of edges')
    plt.tight_layout(pad=0.1)
    if savename:
        plt.savefig(savename,bbox_inches='tight')
    else:
        plt.show()
    plt.close('all')

def run_times_for_cpp(return_run_times_in_dict=False,result_folder='cpp_results'):
    if return_run_times_in_dict:
        d = collections.defaultdict(dict)
    ids = ("arabidopsis",
           "bos",
           "candida",
           "celegans",
           "drosophila",
           "gallus",
           "mus",
           "plasmodium",
           "rattus",
           "sacchcere",
           "sacchpomb"
           )
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    subnet_sizes = []
    number_of_subnets_found = []
    time_fractions = []
    for name in ids:
        for size in sizes:
            filename = result_folder + "/" + name + "_(" + str(size[0]) + "," + str(size[1]) + ").txt"
            try:
                with open(filename,'r') as f:
                    file_is_empty = True
                    for line in f:
                        data = line.split()
                        if data[0] == 'nl-mesu':
                            nl_mesu_time = float(data[1])
                            nl_mesu_number = int(data[2])
                        elif data[0] == 'a-mesu':
                            a_mesu_time = float(data[1])
                            a_mesu_number = int(data[2])
                        file_is_empty = False
                    if file_is_empty:
                        raise Exception # go to except block
                    assert nl_mesu_number == a_mesu_number
                    subnet_sizes.append(str(size))
                    number_of_subnets_found.append(nl_mesu_number)
                    time_fractions.append(a_mesu_time/nl_mesu_time)
                    if return_run_times_in_dict:
                        d[name][size] = (nl_mesu_time,a_mesu_time,nl_mesu_number)
            except:
                if return_run_times_in_dict:
                    d[name][size] = None
                else:
                    pass
    if return_run_times_in_dict:
        return d
    else:
        return subnet_sizes,number_of_subnets_found,time_fractions

def plot_run_times_for_cpp():
    subnet_sizes,number_of_subnets_found,time_fractions = run_times_for_cpp()
    #sc = plt.scatter(subnet_sizes,number_of_subnets_found,s=200,c=time_fractions,norm=colors.TwoSlopeNorm(vcenter=1.0),cmap='BrBG_r',edgecolors='black')
    sc = plt.scatter(subnet_sizes,number_of_subnets_found,s=200,c=time_fractions,cmap='Oranges',edgecolors='black')
    offsets = sc.get_offsets()
    offsets[:,0] += np.random.uniform(-0.1,0.1,offsets.shape[0])
    sc.set_offsets(offsets)
    plt.gca().set_yscale('log')
    cb = plt.colorbar(sc)
    cb.ax.set_ylabel(r'$T_{A-MESU}/T_{NL-MESU}$', rotation=-90, va="bottom")
    plt.xlabel('Subnet size')
    plt.ylabel('Number of subnets found')
    plt.title('C++ relative running times for ppi data')
    plt.tight_layout(pad=0.1)
    plt.savefig('cpp_figures/cpp_relative_run_times.pdf',bbox_inches='tight')

def plot_absolute_running_times_for_cpp(legend=False):
    if legend:
        savename = 'cpp_figures/absolute_running_times.pdf'
    else:
        savename = 'cpp_figures/absolute_running_times_no_legend.pdf'
    ids = ("arabidopsis",
           "bos",
           "candida",
           "celegans",
           "drosophila",
           "gallus",
           "mus",
           "plasmodium",
           "rattus",
           "sacchcere",
           "sacchpomb"
           )
    # change order of sizes
    sizes = ((2,2),(2,3),(3,2),(3,3),(4,2),(4,3))
    colors = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"]
    d = run_times_for_cpp(return_run_times_in_dict=True,result_folder='cpp_results_hammer')
    shared_x_axis = [str(size) for size in sizes]
    individual_y_axes_nl_mesu = []
    individual_y_axes_a_mesu = []
    for name in ids:
        curr_ax_nl = []
        curr_ax_a = []
        for size in sizes:
            curr_ax_nl.append(d[name][size][0] if d[name][size] else None)
            curr_ax_a.append(d[name][size][1] if d[name][size] else None)
        individual_y_axes_nl_mesu.append(curr_ax_nl)
        individual_y_axes_a_mesu.append(curr_ax_a)
    plt.rcParams.update({'lines.linewidth':2})
    fig,ax = plt.subplots(1,1,figsize=(0.79*6.4,0.79*4.8))
    title = ''
    plot_group_of_lines(shared_x_axis, individual_y_axes_nl_mesu, ['-']*len(individual_y_axes_nl_mesu), ids, 'Subnetwork size', r'$t$ (s)', title, plot_to_ax=ax,colors=colors)
    plot_group_of_lines(shared_x_axis, individual_y_axes_a_mesu, ['--']*len(individual_y_axes_a_mesu), ids, 'Subnetwork size', r'$t$ (s)', title, plot_to_ax=ax,colors=colors)
    legend_elements = []
    for ii,name in enumerate(ids):
        legend_elements.append(plt.Line2D([0],[0],marker='o',color=colors[ii],label=name,markerfacecolor=colors[ii],markersize=5,linewidth=0))
    if legend:
        ax.legend(handles=legend_elements,loc='center left',bbox_to_anchor=(1, 0.5))
    if not legend:
        fig.subplots_adjust(top=0.99,right=0.99)
    ax.set_yscale('log')
    ax.set_xlabel('Subnetwork size',size=12)
    ax.set_ylabel(r'$t$ (s)',labelpad=0,size=12)
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    fig.savefig(savename)
    plt.rcParams.update(plt.rcParamsDefault)

def plot_ppi_net_sizes():
    ids = ("arabidopsis",
       "bos",
       "candida",
       "celegans",
       "drosophila",
       "gallus",
       "mus",
       "plasmodium",
       "rattus",
       "sacchcere",
       "sacchpomb"
       )
    n_nodes = [6980,325,367,3879,8215,313,7747,1203,2640,6570,4092]
    n_layers = [7,4,7,6,7,6,7,3,6,7,7]
    n_edges = [19574,360,446,8826,45401,411,21753,2489,4670,304886,69324]
    x_values = [0.02,0.025,-0.025,0,0,0,-0.02,0,0,0,0]
    text_y_multipliers = [0.8,0.75,1.2,1,1,1,1.1,1,1,1,1]
    fig,ax = plt.subplots(figsize=(1.9,4.8))
    colors = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"]
    for ii,ee in enumerate(n_edges):
        ax.scatter(x_values[ii],[ee],c=colors[ii],marker='o',s=50)
        plt.text(0.028,ee*text_y_multipliers[ii],ids[ii],ha='left',size=12)
    ax.set_xlim([-0.04,0.4])
    ax.set_yscale('log')
    ax.set_ylabel('Number of edges',size=12)
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.subplots_adjust(top=0.97,right=1,bottom=0)
    fig.savefig('cpp_figures/edge_numbers.pdf')

def make_table_run_times_for_cpp(filename):
    # table: organism, subnet size, subnet number, nl-mesu time, a-mesu time
    ids = ("arabidopsis",
           "bos",
           "candida",
           "celegans",
           "drosophila",
           "gallus",
           "mus",
           "plasmodium",
           "rattus",
           "sacchcere",
           "sacchpomb"
           )
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    d = run_times_for_cpp(return_run_times_in_dict=True,result_folder='cpp_results_hammer')
    with open(filename,'w') as f:
        for name in ids:
            preamble = name
            for size in sizes:
                subnet_size_str = str(size[0]) + ", " + str(size[1])
                if d[name][size]:
                    number_of_subnets = d[name][size][2]
                    number_of_subnets_str = f"{number_of_subnets:.2e}"
                    nl_mesu_time_str = f"{d[name][size][0]:.2e}"
                    a_mesu_time_str = f"{d[name][size][1]:.2e}"
                    nl_mesu_subnets_per_second_str = f"{number_of_subnets/d[name][size][0]:.2e}" if number_of_subnets else "-"
                    a_mesu_subnets_per_second_str = f"{number_of_subnets/d[name][size][1]:.2e}" if number_of_subnets else "-"
                    # format to LaTeX
                    number_of_subnets_str = "$" + number_of_subnets_str.split('e')[0] + r" \times 10^{" + number_of_subnets_str.split('e')[1].replace('+0','').replace('-0','-') + "}$" if number_of_subnets != 0 else "0"
                    nl_mesu_time_str = "$" + nl_mesu_time_str.split('e')[0] + r" \times 10^{" + nl_mesu_time_str.split('e')[1].replace('+0','').replace('-0','-') + "}$"
                    a_mesu_time_str = "$" + a_mesu_time_str.split('e')[0] + r" \times 10^{" + a_mesu_time_str.split('e')[1].replace('+0','').replace('-0','-') + "}$"
                    nl_mesu_subnets_per_second_str = "$" + nl_mesu_subnets_per_second_str.split('e')[0] + r" \times 10^{" + nl_mesu_subnets_per_second_str.split('e')[1].replace('+0','').replace('-0','-') + "}$" if number_of_subnets else "-"
                    a_mesu_subnets_per_second_str = "$" + a_mesu_subnets_per_second_str.split('e')[0] + r" \times 10^{" + a_mesu_subnets_per_second_str.split('e')[1].replace('+0','').replace('-0','-') + "}$" if number_of_subnets else "-"
                else:
                    number_of_subnets = "-"
                    number_of_subnets_str = "-"
                    nl_mesu_time_str = "-"
                    a_mesu_time_str = "-"
                    nl_mesu_subnets_per_second_str = "-"
                    a_mesu_subnets_per_second_str = "-"
                ending = r" \\" + "\n"
                curr_line = " & ".join([preamble,subnet_size_str,number_of_subnets_str,nl_mesu_time_str,a_mesu_time_str,nl_mesu_subnets_per_second_str,a_mesu_subnets_per_second_str])
                curr_line = curr_line + ending
                #run_times = preamble + " & " + str(size[0]) + ", " + str(size[1]) + " & " + (str(d[name][size][2]) if d[name][size] else "-") + " & " + (str(d[name][size][0]) if d[name][size] else "-") + " & " + (str(d[name][size][1]) if d[name][size] else "-")
                #number_of_subnets = d[name][size][2] if d[name][size] else None
                #subnets_per_second = (str(number_of_subnets/d[name][size][0]) if d[name][size] and number_of_subnets else "-") + " & " + (str(number_of_subnets/d[name][size][1]) if d[name][size] and number_of_subnets else "-")
                #ending = r" \\" + "\n"
                #curr_line = run_times + " & " + subnets_per_second + ending
                f.write(curr_line)
                preamble = ""

##### Model networks benchmarks for cpp

def create_network_in_cpp_format(net_function,**kwargs):
    savename = parse_network_savename(net_function,**kwargs)
    if not os.path.exists(savename):
        M = net_function(**kwargs)
        helpers.save_edgelist_cpp_format(M,savename)
    return savename

def parse_network_savename(net_function,**kwargs):
    savename = 'cpp_benchmark_networks/'
    savename = savename + net_function.__name__
    if kwargs:
        for kw in sorted(kwargs):
            savename = savename + '_' + kw + '=' + str(kwargs[kw]).replace(' ', '')
    return savename

def parse_output_savename(network_inputfilename, subnet_size, save_folder='cpp_benchmark_results/'):
    base_name = network_inputfilename.split('/')[1]
    #save_folder = 'cpp_benchmark_results/'
    return save_folder + base_name + '_' + str(subnet_size).replace(' ','')

def make_er_nets_changing_aspects_generator():
    mean_degree = 3
    n_nodelayers = 1000
    # 1000 nodelayers : p = 0.003 to get <k> = 3
    p = mean_degree/float(n_nodelayers) # approx. Actually should be (n-1) but this is nicer.
    for l in [(40,25),(10,10,10),(8,5,5,5),(2,4,5,5,5)]:
        return_dict = dict()
        kws = dict()
        kws['p'] = p
        kws['l'] = l
        return_dict['net_function'] = helpers.er_multilayer_any_aspects_deg_1_or_greater
        return_dict['kwargs'] = kws
        return_dict['subnet_size'] = (2,)*len(l)
        yield return_dict

def make_geo_mplex_generator():
    # increase nlayers, average degree
    mean_degree = 3
    nnodes = 1000
    nlayers = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30]
    subnet_sizes = [(2,2),(2,3),(3,2),(3,3)]
    for nl in nlayers:
        for subnet_size in subnet_sizes:
            return_dict = dict()
            kws = dict()
            kws['n'] = nnodes
            kws['edges'] = [int((nnodes*mean_degree)/2)]*nl
            kws['couplings'] = 'categorical'
            return_dict['net_function'] = helpers.pn.models.geo
            return_dict['kwargs'] = kws
            return_dict['subnet_size'] = subnet_size
            yield return_dict

def make_geo_mlayer_generator(layers_in_second_aspect=1,nnodes=1000):
    mean_degree_inside = 3
    mean_degree_between = 2
    nlayers_in_first_aspect = [3,4,5,6,7,8,9,10]
    if layers_in_second_aspect > 1:
        subnet_sizes = [(2,2,2),(2,3,2),(3,2,2),(3,3,2)]
    elif layers_in_second_aspect == 1:
        subnet_sizes = [(2,2,1),(2,3,1),(3,2,1),(3,3,1)]
    for nl_first in nlayers_in_first_aspect:
        for subnet_size in subnet_sizes:
            return_dict = dict()
            kws = dict()
            kws['l'] = (nnodes,nl_first,layers_in_second_aspect)
            kws['edges_in_layers'] = int((nnodes*mean_degree_inside)/2)
            kws['edges_between_layers'] = int((nnodes*mean_degree_between)/2)
            return_dict['net_function'] = helpers.geo_multilayer_any_aspects
            return_dict['kwargs'] = kws
            return_dict['subnet_size'] = subnet_size
            yield return_dict

def make_geo_er_mixed_mlayer_generator(nnodes=100,p_er_intra=0.01,p_er_inter=0.01):
    # runs second aspect 1 to 4
    mean_degree_inside = 3
    mean_degree_between = 2
    nlayers_in_first_aspect = [5,10,15,20]
    nlayers_in_second_aspect = [1,2,3,4]
    for layers_in_second_aspect in nlayers_in_second_aspect:
        if layers_in_second_aspect > 1:
            subnet_sizes = [(2,2,2),(2,3,2),(3,2,2),(3,3,2)]
        elif layers_in_second_aspect == 1:
            subnet_sizes = [(2,2,1),(2,3,1),(3,2,1),(3,3,1)]
        for nl_first in nlayers_in_first_aspect:
            for subnet_size in subnet_sizes:
                return_dict = dict()
                kws = dict()
                kws['l'] = (nnodes,nl_first,layers_in_second_aspect)
                kws['edges_in_layers'] = int((nnodes*mean_degree_inside)/2)
                kws['edges_between_layers'] = int((nnodes*mean_degree_between)/2)
                kws['p_er_intra'] = p_er_intra
                kws['p_er_inter'] = p_er_inter
                return_dict['net_function'] = helpers.geo_er_mixed_multilayer_any_aspects
                return_dict['kwargs'] = kws
                return_dict['subnet_size'] = subnet_size
                yield return_dict

def make_geo_er_mixed_mlayer_generator_param_mean_degrees(nnodes,d_geo_intra,d_geo_inter,d_er_intra,d_er_inter):
    nlayers_in_first_aspect = [5,10,15,20]
    #nlayers_in_second_aspect = [1,2]
    nlayers_in_second_aspect = [1,2,3,4]
    for layers_in_second_aspect in nlayers_in_second_aspect:
        if layers_in_second_aspect > 1:
            subnet_sizes = [(2,2,2),(2,3,2),(3,2,2),(3,3,2)]
        elif layers_in_second_aspect == 1:
            subnet_sizes = [(2,2,1),(2,3,1),(3,2,1),(3,3,1)]
        for nl_first in nlayers_in_first_aspect:
            for subnet_size in subnet_sizes:
                return_dict = dict()
                kws = dict()
                kws['l'] = (nnodes,nl_first,layers_in_second_aspect)
                kws['d_geo_intra'] = d_geo_intra
                kws['d_geo_inter'] = d_geo_inter
                kws['d_er_intra'] = d_er_intra
                kws['d_er_inter'] = d_er_inter
                return_dict['net_function'] = helpers.geo_er_mixed_multilayer_any_aspects_mean_deg_parametrization
                return_dict['kwargs'] = kws
                return_dict['subnet_size'] = subnet_size
                yield return_dict

def make_er_mlayer_single_aspect_generator():
    mean_degree = 3
    nlayers = [3,4,5,6,7,8,9,10]
    subnet_sizes = [(2,2),(2,3),(3,2),(3,3)]
    for nl in nlayers:
        for subnet_size in subnet_sizes:
            return_dict = dict()
            kws = dict()
            kws['l'] = (1000,nl)
            # p such that average degree is mean_degree (coupling edges also random)
            # kws['p'] = mean_degree/float(1000*nl)
            # p such that average intralayer degree is mean_degree, and additionally there are other edges
            kws['p'] = mean_degree/1000.0
            return_dict['net_function'] = helpers.er_multilayer_any_aspects_deg_1_or_greater
            return_dict['kwargs'] = kws
            return_dict['subnet_size'] = subnet_size
            yield return_dict

def make_er_mplex_generator():
    # use pymnet with number of edges fixed
    mean_degree = 3
    nnodes = 1000
    nlayers = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30]
    subnet_sizes = [(2,2),(2,3),(3,2),(3,3)]
    for nl in nlayers:
        for subnet_size in subnet_sizes:
            return_dict = dict()
            kws = dict()
            kws['n'] = nnodes
            kws['edges'] = [int((nnodes*mean_degree)/2)]*nl
            return_dict['net_function'] = helpers.pn.models.er
            return_dict['kwargs'] = kws
            return_dict['subnet_size'] = subnet_size
            yield return_dict

def run_benchmark_models_cpp(net_kw_subnet_generator):
    for param_dict in net_kw_subnet_generator:
        inputfile = create_network_in_cpp_format(param_dict['net_function'], **param_dict['kwargs'])
        outputfile = parse_output_savename(inputfile, param_dict['subnet_size'])
        n_aspects = len(param_dict['subnet_size']) - 1
        call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(param_dict['subnet_size']).replace(' ','').strip('()') + "'"
        if not os.path.exists(outputfile):
            os.system(call_str)
            #print(call_str)

def create_networks_without_running_algorithms(net_kw_subnet_generator):
    for param_dict in net_kw_subnet_generator:
        create_network_in_cpp_format(param_dict['net_function'], **param_dict['kwargs'])

##### Extra benchmarks with aggregation method

def run_benchmark_aggregated_cpp(net_kw_subnet_generator):
    # results go to cpp_benchmark_results_aggregated_algo
    for param_dict in net_kw_subnet_generator:
        inputfile = create_network_in_cpp_format(param_dict['net_function'], **param_dict['kwargs'])
        outputfile = parse_output_savename(inputfile, param_dict['subnet_size'], save_folder='cpp_benchmark_results_aggregated_algo/')
        n_aspects = len(param_dict['subnet_size']) - 1
        call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(param_dict['subnet_size']).replace(' ','').strip('()') + "'"
        # add output_method
        call_str = call_str + ' time '
        # add agg algo
        call_str = call_str + ' aggregated'
        if not os.path.exists(outputfile):
            os.system(call_str)

def run_benchmark_aggregated_convenience_script_2aspect_geo_fullnodes():
    net_kw_subnet_generator_list = []
    for nn in (1000,1500,2000,2500,3000,3500,4000,4500,5000,10000):
    #for nn in (1000,):
        for ll in (1,2,3,4):
            # count 3 and 4 layers only for 1000 nodes (very incomplete data for others)
            if ll > 2 and nn > 1000:
                continue
            net_kw_subnet_generator_list.append(make_geo_mlayer_generator(layers_in_second_aspect=ll,nnodes=nn))
    for net_kw_subnet_generator in net_kw_subnet_generator_list:
        run_benchmark_aggregated_cpp(net_kw_subnet_generator)

def run_benchmark_aggregated_convenience_script_2aspect_geo_1000nodes():
    net_kw_subnet_generator_list = []
    #for nn in (1000,1500,2000,2500,3000,3500,4000,4500,5000,10000):
    for nn in (1000,):
        for ll in (1,2,3,4):
            # count 3 and 4 layers only for 1000 nodes (very incomplete data for others)
            if ll > 2 and nn > 1000:
                continue
            net_kw_subnet_generator_list.append(make_geo_mlayer_generator(layers_in_second_aspect=ll,nnodes=nn))
    for net_kw_subnet_generator in net_kw_subnet_generator_list:
        run_benchmark_aggregated_cpp(net_kw_subnet_generator)

def make_shattered_network(nn=100):
    # each layer contains only one connected pair, but aggregated you get the full network
    # also manually insert one 3,2 graphlet into first layer pair
    shattered_net = pn.MultilayerNetwork(aspects=1,fullyInterconnected=False)
    layer = 0
    for node_ii in range(0,nn):
        for node_jj in range(node_ii+1,nn):
            shattered_net[node_ii,layer][node_jj,layer] = 1
            layer = layer + 1
    # add one graphlet
    shattered_net[0,0][0,1] = 1
    savefolder = 'cpp_benchmark_networks_shattered/'
    savename = savefolder+'shattered_'+str(nn)+'_3_2_graphlet.edges'
    helpers.save_edgelist_cpp_format(shattered_net,savename,include_isolated_nls=False)

##### Model network benchmark plotting for cpp

def plot_group_of_lines(shared_x_axis,individual_y_axes,formats,line_labels,x_axis_label,y_axis_label,title,plot_to_ax=None,colors=None):
    if plot_to_ax is None:
        fig,ax = plt.subplots()
    else:
        ax = plot_to_ax
        fig = ax.get_figure()
    for ii in range(len(individual_y_axes)):
        if colors:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii],color=colors[ii])
        else:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii])
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(y_axis_label)
    ax.set_title(title)
    ax.set_xticks(shared_x_axis)
    return fig,ax

def plot_group_of_lines_xlabels(shared_x_axis,x_axis_labels,individual_y_axes,formats,line_labels,x_axis_label,y_axis_label,title,plot_to_ax=None,colors=None):
    if plot_to_ax is None:
        fig,ax = plt.subplots()
    else:
        ax = plot_to_ax
        fig = ax.get_figure()
    for ii in range(len(individual_y_axes)):
        if colors:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii],color=colors[ii])
        else:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii])
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(y_axis_label)
    ax.set_title(title)
    ax.set_xticks(ticks=shared_x_axis,labels=x_axis_labels)
    return fig,ax

def plot_mplex_relative_vs_net_size(net_kw_subnet_generator=make_geo_mplex_generator(),savename='cpp_benchmark_figures/geo_mplex.pdf',title='GEO mplex',plot_to_ax=None):
    #shared_x_axis = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    shared_x_axis = [3,4,5,6,7,8,9,10]
    res_dict = dict()
    for param_dict in net_kw_subnet_generator:
        try:
            nl = len(param_dict['kwargs']['edges'])
        except:
            nl = param_dict['kwargs']['l'][1] # number of layers in first aspect
        outputfile = parse_output_savename(parse_network_savename(param_dict['net_function'],**param_dict['kwargs']),param_dict['subnet_size'])
        try:
            nl_mesu_res,a_mesu_res = read_temp_file_res(outputfile)
        except:
            nl_mesu_res,a_mesu_res = ((None,None),(None,None))
        assert(nl_mesu_res[1] == a_mesu_res[1])
        # default value is list of -1's with length len(shared_x_axis). Insert into correct place according to nl
        # value: t(a-mesu)/t(nl-mesu)
        try:
            res_dict.setdefault(param_dict['subnet_size'],[-1]*len(shared_x_axis))[shared_x_axis.index(nl)] = a_mesu_res[0]/nl_mesu_res[0]
        except:
            res_dict.setdefault(param_dict['subnet_size'],[-1]*len(shared_x_axis))[shared_x_axis.index(nl)] = None
    individual_y_axes = []
    formats = []
    line_labels = []
    for kk in sorted(res_dict.keys()):
        individual_y_axes.append(res_dict[kk])
        formats.append('-')
        line_labels.append(str(kk))
    fig,ax = plot_group_of_lines(shared_x_axis,individual_y_axes,formats,line_labels,'Number of layers',r'$t_{a-mesu}/t_{nl-mesu}$',title,plot_to_ax=plot_to_ax)
    ax.set_yscale('log')
    ax.legend()
    ax.plot(shared_x_axis,[1]*len(shared_x_axis),'--k')
    if not plot_to_ax:
        fig.savefig(savename)

def plot_absolute_times_vs_net_size(net_kw_subnet_generator,savename,title):
    #shared_x_axis = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    shared_x_axis = [3,4,5,6,7,8,9,10]
    res_dict_with_separate_a_and_nl = dict()
    res_dict_with_separate_a_and_nl['a-mesu'] = dict()
    res_dict_with_separate_a_and_nl['nl-mesu'] = dict()
    for param_dict in net_kw_subnet_generator:
        try:
            nl = len(param_dict['kwargs']['edges'])
        except:
            nl = param_dict['kwargs']['l'][1] # number of layers in first aspect
        outputfile = parse_output_savename(parse_network_savename(param_dict['net_function'],**param_dict['kwargs']),param_dict['subnet_size'])
        try:
            nl_mesu_res,a_mesu_res = read_temp_file_res(outputfile)
        except:
            nl_mesu_res,a_mesu_res = ((None,None),(None,None))
        assert(nl_mesu_res[1] == a_mesu_res[1])
        # default value is list of -1's with length len(shared_x_axis). Insert into correct place according to nl
        # value: t(a-mesu)/t(nl-mesu)
        try:
            res_dict_with_separate_a_and_nl['a-mesu'].setdefault(param_dict['subnet_size'],[-1]*len(shared_x_axis))[shared_x_axis.index(nl)] = a_mesu_res[0]
        except:
            res_dict_with_separate_a_and_nl['a-mesu'].setdefault(param_dict['subnet_size'],[-1]*len(shared_x_axis))[shared_x_axis.index(nl)] = None
        try:
            res_dict_with_separate_a_and_nl['nl-mesu'].setdefault(param_dict['subnet_size'],[-1]*len(shared_x_axis))[shared_x_axis.index(nl)] = nl_mesu_res[0]
        except:
            res_dict_with_separate_a_and_nl['nl-mesu'].setdefault(param_dict['subnet_size'],[-1]*len(shared_x_axis))[shared_x_axis.index(nl)] = None
    individual_y_axes = []
    formats = []
    line_labels = []
    colors = []
    # default plt color cycle
    color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    for ii,kk in enumerate(sorted(res_dict_with_separate_a_and_nl['nl-mesu'].keys())):
        individual_y_axes.append(res_dict_with_separate_a_and_nl['nl-mesu'][kk])
        formats.append('-')
        line_labels.append(str(kk))
        colors.append(color_cycle[ii])
    for ii,kk in enumerate(sorted(res_dict_with_separate_a_and_nl['a-mesu'].keys())):
        individual_y_axes.append(res_dict_with_separate_a_and_nl['a-mesu'][kk])
        formats.append('--')
        line_labels.append(str(kk))
        colors.append(color_cycle[ii])
    fig,ax = plot_group_of_lines(shared_x_axis,individual_y_axes,formats,line_labels,'Number of layers',r'$t$',title,colors=colors)
    ax.set_yscale('log')
    ax.legend()
    #ax.plot(shared_x_axis,[1]*len(shared_x_axis),'--k')
    fig.savefig(savename)

def plot_scatter_times_vs_number_of_subnets(net_kw_subnet_generator_list,savename):
    nl_mesu_scatter = [[],[]]
    a_mesu_scatter = [[],[]]
    success_list = []
    for net_kw_subnet_generator in net_kw_subnet_generator_list:
        for param_dict in net_kw_subnet_generator:
            outputfile = parse_output_savename(parse_network_savename(param_dict['net_function'],**param_dict['kwargs']),param_dict['subnet_size'])
            try:
                nl_mesu_res,a_mesu_res = read_temp_file_res(outputfile)
                assert(nl_mesu_res[1] == a_mesu_res[1])
                nl_mesu_scatter[0].append(nl_mesu_res[1])
                nl_mesu_scatter[1].append(nl_mesu_res[0])
                a_mesu_scatter[0].append(a_mesu_res[1])
                a_mesu_scatter[1].append(a_mesu_res[0])
                success_list.append(param_dict)
            except:
                pass
    fig,ax = plt.subplots(figsize=(0.7*6.4,0.7*4.8))
    ax.scatter(nl_mesu_scatter[0],nl_mesu_scatter[1],color='darkred',marker='x',label='nlse',linewidth=1)
    ax.scatter(a_mesu_scatter[0],a_mesu_scatter[1],color='darkgreen',marker='+',label='elsse',linewidth=1)
    print(a_mesu_scatter[0] == nl_mesu_scatter[0])
    print(a_mesu_scatter[1] == nl_mesu_scatter[1])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Number of subnetworks',size=12)
    ax.set_ylabel(r'$t$ (s)',size=12)
    ax.legend()
    fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
    if savename:
        fig.savefig(savename)
        return success_list
    else:
        return success_list,fig,ax

def plot_scatter_convenience_script(plot_fit_lines=True):
    net_kw_subnet_generator_list = []
    for nn in (1000,1500,2000,2500,3000,3500,4000,4500,5000,10000):
        for ll in (1,2,3,4):
            # count 3 and 4 layers only for 1000 nodes (very incomplete data for others)
            if ll > 2 and nn > 1000:
                continue
            net_kw_subnet_generator_list.append(make_geo_mlayer_generator(layers_in_second_aspect=ll,nnodes=nn))
    savename = 'cpp_benchmark_figures/absolute_scatter_geo.pdf'
    if plot_fit_lines:
        success_list,fig,ax = plot_scatter_times_vs_number_of_subnets(net_kw_subnet_generator_list, savename=None)
        fit_x_1 = np.logspace(3,6.3,num=50,base=10)
        fit_x_2 = np.logspace(3.40,6.47,num=50,base=10)
        fit_y_1 = []
        fit_y_2 =  []
        for x in fit_x_1:
            fit_y_1.append(10**-1.5 * (x**1.2))
        for x in fit_x_2:
            fit_y_2.append(10**-5 * (x**1))
        ax.plot(fit_x_1,fit_y_1,color='k',linestyle='--')
        plt.annotate('',xy=(1500,300),xytext=(1200,16000),arrowprops=dict(arrowstyle='simple,tail_width=0.1',facecolor='black'))
        plt.text(800,20000,r'$t = 10^{-1.5} \times n_{subnets}^{1.2}$',fontsize=12)
        ax.plot(fit_x_2,fit_y_2,color='k',linestyle='--')
        plt.annotate('',xy=(1000000,6),xytext=(800000,0.25),arrowprops=dict(arrowstyle='simple,tail_width=0.1',facecolor='black'))
        plt.text(120000,0.1,r'$t = 10^{-5} \times n_{subnets}$',fontsize=12)
        fig.savefig(savename)
    else:
        success_list = plot_scatter_times_vs_number_of_subnets(net_kw_subnet_generator_list, savename)
    for kwd in success_list:
        (nnodes,nl_first,layers_in_second_aspect) = kwd['kwargs']['l']
        print(nnodes,nl_first,layers_in_second_aspect)

def plot_combined_convenience_script(case):
    # case 1 : two-aspect multilayer geometric (4 imgs)
    if case == 1:
        #plt.rcParams.update({'legend.fontsize': 8.4,'legend.handlelength': 1,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        plt.rcParams.update({'legend.fontsize': 6,'legend.handlelength': 0.8,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        #fig,axs = plt.subplots(2, 2, sharex='all', sharey='all')
        fig,axs = plt.subplots(2, 2, sharex='all', sharey='all', figsize=(0.7*6.4,0.7*4.8))
        # 1000 n 1 l
        net_kw_subnet_generator = make_geo_mlayer_generator(layers_in_second_aspect=1,nnodes=1000)
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs[0][0])
        axs[0][0].set_xlabel(None)
        axs[0][0].set_ylabel(None)
        axs[0][0].text(3,0.08,'a)',fontsize=12)
        axs[0][0].legend(ncols=4)
        axs[0][0].set_ylim([0.005,9])
        # 1000 n 2 l
        net_kw_subnet_generator = make_geo_mlayer_generator(layers_in_second_aspect=2,nnodes=1000)
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs[0][1])
        axs[0][1].set_xlabel(None)
        axs[0][1].set_ylabel(None)
        axs[0][1].text(3,0.08,'b)',fontsize=12)
        axs[0][1].legend(ncols=4)
        # 10 000 n 1 l
        net_kw_subnet_generator = make_geo_mlayer_generator(layers_in_second_aspect=1,nnodes=10000)
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs[1][0])
        axs[1][0].set_xlabel(None)
        axs[1][0].set_ylabel(None)
        axs[1][0].text(3,0.08,'c)',fontsize=12)
        axs[1][0].legend(ncols=4)
        # 10 000 n 2 l
        net_kw_subnet_generator = make_geo_mlayer_generator(layers_in_second_aspect=2,nnodes=10000)
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs[1][1])
        axs[1][1].set_xlabel(None)
        axs[1][1].set_ylabel(None)
        axs[1][1].text(3,0.08,'d)',fontsize=12)
        axs[1][1].legend(ncols=4)
        fig.subplots_adjust(wspace=0, hspace=0)
        #fig.subplots_adjust(top=0.99,right=0.99)
        fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
        fig.supxlabel('            Number of elementary layers in first aspect')
        fig.supylabel(r'         $t_{elsse}$ / $t_{nlse}$',x=0)
        fig.savefig('cpp_benchmark_figures/geo_mlayer_combined_1000_10000.pdf')
        plt.rcParams.update(plt.rcParamsDefault)
    # case 2 : mplex geo and er
    # set xaxis in plotter to 3...20 !!!
    if case == 2:
        #plt.rcParams.update({'legend.fontsize': 8.4,'legend.handlelength': 1,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        plt.rcParams.update({'legend.fontsize': 6,'legend.handlelength': 0.8,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        #fig,axs = plt.subplots(1, 2, sharex='all', sharey='all',figsize=(6.4,2.4))
        fig,axs = plt.subplots(1, 2, sharex='all', sharey='all',figsize=(0.7*6.4,1.2*0.7*2.4))
        # geo
        net_kw_subnet_generator = make_geo_mplex_generator()
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs[0])
        axs[0].set_xlabel(None)
        axs[0].set_ylabel(None)
        axs[0].text(3,0.4,'a)',fontsize=12)
        axs[0].legend(ncols=4)
        axs[0].set_ylim(0.09,4.2)
        xticklabs = axs[0].get_xticklabels()
        plt.setp(axs[0].get_xticklabels(), visible=False)
        plt.setp(xticklabs[::2], visible=True)
        # er
        net_kw_subnet_generator = make_er_mplex_generator()
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs[1])
        axs[1].set_xlabel(None)
        axs[1].set_ylabel(None)
        axs[1].text(3,0.4,'b)',fontsize=12)
        axs[1].legend(ncols=4)
        xticklabs = axs[1].get_xticklabels()
        plt.setp(axs[1].get_xticklabels(), visible=False)
        plt.setp(xticklabs[::2], visible=True)
        fig.subplots_adjust(wspace=0, hspace=0)
        #fig.subplots_adjust(top=0.99,right=0.99,bottom=0.2)
        fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.22)
        fig.supxlabel('           Number of layers')
        fig.supylabel(r'         $t_{elsse}$ / $t_{nlse}$',x=0)
        fig.savefig('cpp_benchmark_figures/mplex_geo_er_combined.pdf')
        plt.rcParams.update(plt.rcParamsDefault)
    # case 3: supplementary, geo mlayer with 3 in 2nd aspect
    if case == 3:
        #plt.rcParams.update({'legend.fontsize': 8.4,'legend.handlelength': 1,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        plt.rcParams.update({'legend.fontsize': 12,'legend.handlelength': 0.8,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        #fig,axs = plt.subplots(2, 2, sharex='all', sharey='all')
        fig,axs = plt.subplots(1, 1, sharex='all', sharey='all', figsize=(0.7*6.4,0.7*4.8))
        # 1000 n 1 l
        net_kw_subnet_generator = make_geo_mlayer_generator(layers_in_second_aspect=3,nnodes=1000)
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs)
        axs.set_xlabel(None)
        axs.set_ylabel(None)
        axs.text(3,0.08,'a)',fontsize=12)
        axs.legend(ncols=4)
        axs.set_ylim([0.005,4])
        #fig.subplots_adjust(top=0.99,right=0.99)
        fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
        fig.supxlabel('            Number of elementary layers in first aspect')
        fig.supylabel(r'         $t_{elsse}$ / $t_{nlse}$',x=0)
        fig.savefig('cpp_benchmark_figures/geo_mlayer_3_1000n.pdf')
        plt.rcParams.update(plt.rcParamsDefault)
    # case 4: supplementary, geo mlayer with 4 in 2nd aspect
    if case == 4:
        #plt.rcParams.update({'legend.fontsize': 8.4,'legend.handlelength': 1,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        plt.rcParams.update({'legend.fontsize': 12,'legend.handlelength': 0.8,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
        #fig,axs = plt.subplots(2, 2, sharex='all', sharey='all')
        fig,axs = plt.subplots(1, 1, sharex='all', sharey='all', figsize=(0.7*6.4,0.7*4.8))
        # 1000 n 1 l
        net_kw_subnet_generator = make_geo_mlayer_generator(layers_in_second_aspect=4,nnodes=1000)
        plot_mplex_relative_vs_net_size(net_kw_subnet_generator,savename='should_not_exist.pdf',title=None,plot_to_ax=axs)
        axs.set_xlabel(None)
        axs.set_ylabel(None)
        axs.text(3,0.08,'b)',fontsize=12)
        axs.legend(ncols=4)
        axs.set_ylim([0.005,4])
        #fig.subplots_adjust(top=0.99,right=0.99)
        fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
        fig.supxlabel('            Number of elementary layers in first aspect')
        fig.supylabel(r'         $t_{elsse}$ / $t_{nlse}$',x=0)
        fig.savefig('cpp_benchmark_figures/geo_mlayer_4_1000n.pdf')
        plt.rcParams.update(plt.rcParamsDefault)

def load_benchmark_result_from_all_algos(param_dict,return_outputfiles=False):
    outputfiles = []
    # nl-mesu and a-mesu
    outputfiles.append(parse_output_savename(parse_network_savename(param_dict['net_function'],**param_dict['kwargs']),param_dict['subnet_size'],save_folder='cpp_benchmark_results/'))
    # aggregated
    outputfiles.append(parse_output_savename(parse_network_savename(param_dict['net_function'],**param_dict['kwargs']),param_dict['subnet_size'],save_folder='cpp_benchmark_results_aggregated_algo/'))
    # defaults to None
    #nl_mesu_time,nl_mesu_number,a_mesu_time,a_mesu_number,agg_time,agg_number = None,None,None,None,None,None
    res_dict = {}
    for key in ['nl-mesu_time','nl-mesu_number','a-mesu_time','a-mesu_number','aggregated_time','aggregated_number']:
        res_dict[key] = None
    for outputfile in outputfiles:
        try:
            with open(outputfile,'r') as f:
                file_is_empty = True
                for line in f:
                    data = line.split()
                    res_dict[data[0]+'_time'] = float(data[1])
                    res_dict[data[0]+'_number'] = int(data[2])
                    file_is_empty = False
                if file_is_empty:
                    raise Exception # go to except block
        except:
            pass
    if return_outputfiles:
        return res_dict,outputfiles
    else:
        return res_dict


def nested_dict():
    return collections.defaultdict(nested_dict)

def load_benchmark_result_by_subnet_size_all_algos(net_kw_subnet_generator,return_outputfiles=False):
    res_dict = nested_dict()
    if return_outputfiles:
        outputfile_dict = nested_dict()
    for param_dict in net_kw_subnet_generator:
        nnodes,nlayers_in_first_aspect,nlayers_in_second_aspect = param_dict['kwargs']['l']
        res,outputfiles =  load_benchmark_result_from_all_algos(param_dict,return_outputfiles=True)
        res_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect][param_dict['subnet_size']] = res
        if return_outputfiles:
            outputfile_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect][param_dict['subnet_size']] = outputfiles
    if return_outputfiles:
        return res_dict,outputfile_dict
    else:
        return res_dict

def load_benchmark_result_from_mplex_all_algos(net_kw_subnet_generator):
    res_dict = nested_dict()
    for param_dict in net_kw_subnet_generator:
        nnodes = param_dict['kwargs']['n']
        nlayers = len(param_dict['kwargs']['edges'])
        res_dict[nnodes][nlayers][param_dict['subnet_size']] = load_benchmark_result_from_all_algos(param_dict)
    return res_dict

def plot_geo_er_mixed_all_algos_separate_plot(full_dict,title,nlayers_in_first_aspect=5,nlayers_in_second_aspect=1,nnodes=100):
    formats = {'nl-mesu' : '-', 'a-mesu' : '--', 'aggregated' : '-.'}
    #full_dict = load_benchmark_result_by_subnet_size_all_algos(net_kw_subnet_generator)
    time_dict = full_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect]
    x_axis_labels_tuples = sorted(time_dict.keys(),key=lambda x : x[-1])
    x_axis_labels = [str(x) for x in x_axis_labels_tuples]
    shared_x_axis = range(len(x_axis_labels))
    # nl-mesu, a-mesu, aggregated
    all_y_vals = [[],[],[]]
    all_formats = [formats['nl-mesu'],formats['a-mesu'],formats['aggregated']]
    all_labels = ['nlse','elsse','ad-esu']
    all_colors = ['#000000','#000000','#000000']
    for subnet_size in x_axis_labels_tuples:
        all_y_vals[0].append(time_dict[subnet_size]['nl-mesu_time'])
        all_y_vals[1].append(time_dict[subnet_size]['a-mesu_time'])
        all_y_vals[2].append(time_dict[subnet_size]['aggregated_time'])
    # dirty hack for avoiding pdf broken errors with all None values
    if any(all_y_vals[0]) or any(all_y_vals[1]) or any(all_y_vals[2]):
        fig,ax = plot_group_of_lines_xlabels(shared_x_axis, x_axis_labels, individual_y_axes=all_y_vals, formats=all_formats, line_labels=all_labels, x_axis_label='Subnetwork size', y_axis_label='t (s)', title=title, colors=all_colors)
        ax.set_yscale('log')
        ax.legend(loc='upper left')
        return fig,ax
    else:
        return None,None

def plot_all_individual_geo_er_convenience_script():
    nn = 100
    savefolder = 'cpp_benchmark_figures/geo_er_mixed_'+str(nn)+'nodes/'
    nl1 = [5,10,15,20]
    nl2 = [1,2,3,4]
    equal_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,1.5,1.5,1.5,1.5)
    geo_1_er_2_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,d_geo_intra=1,d_geo_inter=1,d_er_intra=2,d_er_inter=2)
    geo_2_er_1_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,d_geo_intra=2,d_geo_inter=2,d_er_intra=1,d_er_inter=1)
    gens = [equal_generator,geo_1_er_2_generator,geo_2_er_1_generator]
    full_dicts = []
    for gen in gens:
        full_dicts.append(load_benchmark_result_by_subnet_size_all_algos(gen))
    gen_names = ['geo1.5_er1.5','geo1_er2','geo2_er1']
    for full_dict,gname in zip(full_dicts,gen_names):
        for nlayers_in_first_aspect in nl1:
            for nlayers_in_second_aspect in nl2:
                title = '1.asp: '+str(nlayers_in_first_aspect)+' 2.asp: '+str(nlayers_in_second_aspect)+' '+gname
                #try:
                fig,ax = plot_geo_er_mixed_all_algos_separate_plot(full_dict,title,nlayers_in_first_aspect,nlayers_in_second_aspect,nn)
                if fig is not None:
                    fig.savefig(savefolder+title+'.pdf')
                #except:
                #    pass
    plt.close('all')

def plot_geo_er_mixed_absolute_vs_subnet_size(net_kw_subnet_generator,nlayers_in_first_aspect_list,nlayers_in_second_aspect,plot_to_ax):
    #color_dict = {5:'#c2e699',10:'#78c679',15:'#31a354',20:'#006837'}
    color_dict = {5:'#D81B60',10:'#004D40',15:'#DEAB0F',20:'#1E88E5'}
    formats = {'nl-mesu' : '-', 'a-mesu' : '--', 'aggregated' : ':'}
    nnodes = 100
    #nl1 = [5,10,15,20]
    #nl1 = [10,20]
    nl1 = nlayers_in_first_aspect_list
    full_dict = load_benchmark_result_by_subnet_size_all_algos(net_kw_subnet_generator)
    time_dict = full_dict[nnodes][nlayers_in_second_aspect]
    if nlayers_in_second_aspect == 1:
        subnet_sizes = [(2,2,1),(2,3,1),(3,2,1),(3,3,1)]
    else:
        subnet_sizes = [(2,2,2),(2,3,2),(3,2,2),(3,3,2)]
    x_axis_labels = [''.join(str(x).split()) for x in subnet_sizes]
    shared_x_axis = range(len(x_axis_labels))
    #fig,ax = plt.subplots(1,1)
    for nlayers_in_first_aspect in nl1:
        nl1_y_vals = [[],[],[]]
        for subnet_size in subnet_sizes:
            nl1_y_vals[0].append(float(time_dict[nlayers_in_first_aspect][subnet_size]['nl-mesu_time']))
            nl1_y_vals[1].append(float(time_dict[nlayers_in_first_aspect][subnet_size]['a-mesu_time']))
            nl1_y_vals[2].append(float(time_dict[nlayers_in_first_aspect][subnet_size]['aggregated_time']))
        nl1_formats = [formats['nl-mesu'],formats['a-mesu'],formats['aggregated']]
        nl1_colors = [color_dict[nlayers_in_first_aspect],color_dict[nlayers_in_first_aspect],color_dict[nlayers_in_first_aspect]]
        nl1_labels = ['nlse '+str(nlayers_in_first_aspect),'elsse '+str(nlayers_in_first_aspect),'ad-esu '+str(nlayers_in_first_aspect)]
        plot_group_of_lines_xlabels(shared_x_axis, x_axis_labels, individual_y_axes=nl1_y_vals, formats=nl1_formats, line_labels=nl1_labels, x_axis_label='Subnetwork size', y_axis_label='t (s)', title=None, colors=nl1_colors, plot_to_ax=plot_to_ax)
    plot_to_ax.set_yscale('log')

def plot_geo_er_mixed_convenience_script():
    color_dict = {5:'#D81B60',10:'#004D40',15:'#DEAB0F',20:'#1E88E5'}
    nl1 = [5,20]
    nl1 = [5,10,15,20]
    #plt.rcParams.update({'legend.fontsize': 8.4,'legend.handlelength': 1,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
    #plt.rcParams.update({'legend.fontsize': 6,'legend.handlelength': 0.8,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
    plt.rcParams.update({'legend.fontsize': 8,'legend.handlelength': 3.5,'legend.loc':'lower left','legend.columnspacing': 0.7,'legend.handletextpad': 0.3,'lines.linewidth':2})
    #fig,axs = plt.subplots(2, 2, sharex='all', sharey='all')
    #fig,axs = plt.subplots(2, 2, sharex='col', sharey='all', figsize=(0.7*6.4,0.7*4.8))
    fig,axs = plt.subplots(2, 2, sharex='col', sharey='all', figsize=(0.7*6.4,0.7*6))
    # geo 2 er 1
    # 1 l in 2nd aspect
    net_kw_subnet_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,2,2,1,1)
    plot_geo_er_mixed_absolute_vs_subnet_size(net_kw_subnet_generator,nl1,1,plot_to_ax=axs[0][0])
    axs[0][0].set_xlabel(None)
    axs[0][0].set_ylabel(None)
    axs[0][0].text(0,2000,'a)',fontsize=12)
    #axs[0][0].legend(ncols=4)
    axs[0][0].set_ylim([0.005,20000])
    # geo 2 er 1
    # 4 l in second aspect
    net_kw_subnet_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,2,2,1,1)
    plot_geo_er_mixed_absolute_vs_subnet_size(net_kw_subnet_generator,nl1,4,plot_to_ax=axs[0][1])
    axs[0][1].set_xlabel(None)
    axs[0][1].set_ylabel(None)
    axs[0][1].text(0,2000,'b)',fontsize=12)
    #axs[0][1].legend(ncols=4)
    # geo 1 er 2
    # 1 l in 2nd aspect
    net_kw_subnet_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,1,1,2,2)
    plot_geo_er_mixed_absolute_vs_subnet_size(net_kw_subnet_generator,nl1,1,plot_to_ax=axs[1][0])
    axs[1][0].set_xlabel(None)
    axs[1][0].set_ylabel(None)
    axs[1][0].text(0,2000,'c)',fontsize=12)
    #axs[1][0].legend(ncols=4)
    # geo 1 er 2
    # 4 l in 2nd aspect
    net_kw_subnet_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,1,1,2,2)
    plot_geo_er_mixed_absolute_vs_subnet_size(net_kw_subnet_generator,nl1,4,plot_to_ax=axs[1][1])
    axs[1][1].set_xlabel(None)
    axs[1][1].set_ylabel(None)
    axs[1][1].text(0,2000,'d)',fontsize=12)
    #axs[1][1].legend(ncols=4)
    #axs[0][0].tick_params(axis='x', labelrotation=45)
    #axs[0][1].tick_params(axis='x', labelrotation=45)
    axs[1][0].tick_params(axis='x', labelrotation=75)
    axs[1][1].tick_params(axis='x', labelrotation=75)
    fig.subplots_adjust(wspace=0, hspace=0)
    #fig.subplots_adjust(top=0.99,right=0.99)
    #fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
    fig.subplots_adjust(top=0.95,right=0.99,left=0.13,bottom=0.18)
    fig.supxlabel('          Subnetwork size',y=0)
    fig.supylabel(r'           $t$ (s)',x=0)
    #fig.legend(*axs[0][1].get_legend_handles_labels(),loc='lower right',bbox_to_anchor=(2, 0),fancybox=False,shadow=False,ncol=4)
    lines2 = [Line2D([0], [1], color='k', linewidth=2, linestyle=ls) for ls in ['-','--',':']]
    labels2 = ['nlse','elsse','ad-esu']
    legend2 = plt.legend(lines2, labels2, loc=8, ncols=3, bbox_to_anchor=(0,2),frameon=False,fancybox=False,shadow=False)
    if len(nl1) == 2:
        elem_layer_text_start_point = 0.25
        under10_addition = 0.03
        over10_addition = 0.05
    else:
        elem_layer_text_start_point = 0.196
        under10_addition = 0.032
        over10_addition = 0.05
    addition = 0
    fig.text(elem_layer_text_start_point,0.955,r'Number of elementary layers in first aspect:',fontsize=8)
    for nlayers_in_first_aspect in nl1:
        fig.text(elem_layer_text_start_point+0.56+addition,0.955,str(nlayers_in_first_aspect),color=color_dict[nlayers_in_first_aspect],fontsize=8,weight='bold')
        if nlayers_in_first_aspect < 10:
            addition = addition + under10_addition
        else:
            addition = addition + over10_addition
    if len(nl1) == 2:
        fig.savefig('cpp_benchmark_figures/geo_er_mixed_combined_2_1_1_2_1asp_5_20.pdf')
    else:
        fig.savefig('cpp_benchmark_figures/geo_er_mixed_combined_2_1_1_2_1asp_all.pdf')
    plt.rcParams.update(plt.rcParamsDefault)

def plot_mplex_geo_and_er_absolute_vs_subnet_size(net_kw_subnet_generator,nlayers_list,plot_to_ax):
    #color_dict = {5:'#c2e699',10:'#78c679',15:'#31a354',20:'#006837'}
    color_dict = {5:'#D81B60',10:'#004D40',15:'#DEAB0F',20:'#1E88E5',30:'#D55E00'}
    formats = {'nl-mesu' : '-', 'a-mesu' : '--', 'aggregated' : ':'}
    nnodes = 1000
    #nl1 = [5,10,15,20]
    #nl1 = [10,20]
    nl1 = nlayers_list
    full_dict = load_benchmark_result_from_mplex_all_algos(net_kw_subnet_generator)
    time_dict = full_dict[nnodes]
    subnet_sizes = [(2,2),(2,3),(3,2),(3,3)]
    x_axis_labels = [''.join(str(x).split()) for x in subnet_sizes]
    shared_x_axis = range(len(x_axis_labels))
    #fig,ax = plt.subplots(1,1)
    for nlayers_in_first_aspect in nl1:
        nl1_y_vals = [[],[],[]]
        for subnet_size in subnet_sizes:
            nl1_y_vals[0].append(float(time_dict[nlayers_in_first_aspect][subnet_size]['nl-mesu_time']))
            nl1_y_vals[1].append(float(time_dict[nlayers_in_first_aspect][subnet_size]['a-mesu_time']))
            nl1_y_vals[2].append(float(time_dict[nlayers_in_first_aspect][subnet_size]['aggregated_time']))
        nl1_formats = [formats['nl-mesu'],formats['a-mesu'],formats['aggregated']]
        nl1_colors = [color_dict[nlayers_in_first_aspect],color_dict[nlayers_in_first_aspect],color_dict[nlayers_in_first_aspect]]
        nl1_labels = ['nlse '+str(nlayers_in_first_aspect),'elsse '+str(nlayers_in_first_aspect),'ad-esu '+str(nlayers_in_first_aspect)]
        plot_group_of_lines_xlabels(shared_x_axis, x_axis_labels, individual_y_axes=nl1_y_vals, formats=nl1_formats, line_labels=nl1_labels, x_axis_label='Subnetwork size', y_axis_label='t (s)', title=None, colors=nl1_colors, plot_to_ax=plot_to_ax)
    plot_to_ax.set_yscale('log')

def plot_mplex_geo_er_30l_convenience_script():
    color_dict = {5:'#D81B60',10:'#004D40',15:'#DEAB0F',20:'#1E88E5',30:'#D55E00'}
    nl1 = [10,20,30]
    #plt.rcParams.update({'legend.fontsize': 8.4,'legend.handlelength': 1,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
    #plt.rcParams.update({'legend.fontsize': 6,'legend.handlelength': 0.8,'legend.loc':'lower left','legend.columnspacing': 0.4,'legend.handletextpad': 0.2,'lines.linewidth':2})
    plt.rcParams.update({'legend.fontsize': 8,'legend.handlelength': 3.5,'legend.loc':'lower left','legend.columnspacing': 0.7,'legend.handletextpad': 0.3,'lines.linewidth':2})
    #fig,axs = plt.subplots(2, 2, sharex='all', sharey='all')
    #fig,axs = plt.subplots(2, 2, sharex='col', sharey='all', figsize=(0.7*6.4,0.7*4.8))
    fig,axs = plt.subplots(1, 2, sharex='col', sharey='all', figsize=(0.7*6.4,0.7*6*0.6))
    # geo
    net_kw_subnet_generator = make_geo_mplex_generator()
    plot_mplex_geo_and_er_absolute_vs_subnet_size(net_kw_subnet_generator,nl1,plot_to_ax=axs[0])
    axs[0].set_xlabel(None)
    axs[0].set_ylabel(None)
    axs[0].text(0,2000,'a)',fontsize=12)
    #axs[0][0].legend(ncols=4)
    #axs[0].set_ylim([0.005,20000])
    # er
    net_kw_subnet_generator = make_er_mplex_generator()
    plot_mplex_geo_and_er_absolute_vs_subnet_size(net_kw_subnet_generator,nl1,plot_to_ax=axs[1])
    axs[1].set_xlabel(None)
    axs[1].set_ylabel(None)
    axs[1].text(0,2000,'b)',fontsize=12)
    fig.subplots_adjust(wspace=0, hspace=0)
    #fig.subplots_adjust(top=0.99,right=0.99)
    #fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
    fig.subplots_adjust(top=0.918,right=0.99,left=0.13,bottom=0.25)
    fig.supxlabel('          Subnetwork size',y=0)
    fig.supylabel(r'           $t$ (s)',x=0)
    axs[0].tick_params(axis='x', labelrotation=75)
    axs[1].tick_params(axis='x', labelrotation=75)
    #fig.legend(*axs[0][1].get_legend_handles_labels(),loc='lower right',bbox_to_anchor=(2, 0),fancybox=False,shadow=False,ncol=4)
    lines2 = [Line2D([0], [1], color='k', linewidth=2, linestyle=ls) for ls in ['-','--',':']]
    labels2 = ['nlse','elsse','ad-esu']
    legend2 = plt.legend(lines2, labels2, loc=8, ncols=3, bbox_to_anchor=(0,0.999),frameon=False,fancybox=False,shadow=False)
    if len(nl1) == 2:
        elem_layer_text_start_point = 0.25
        under10_addition = 0.03
        over10_addition = 0.05
    else:
        elem_layer_text_start_point = 0.196
        under10_addition = 0.032
        over10_addition = 0.05
    addition = 0
    fig.text(elem_layer_text_start_point,0.926,r'Number of elementary layers in first aspect:',fontsize=8)
    for nlayers_in_first_aspect in nl1:
        fig.text(elem_layer_text_start_point+0.56+addition,0.926,str(nlayers_in_first_aspect),color=color_dict[nlayers_in_first_aspect],fontsize=8,weight='bold')
        if nlayers_in_first_aspect < 10:
            addition = addition + under10_addition
        else:
            addition = addition + over10_addition
    fig.savefig('cpp_benchmark_figures/mplex_geo_er_absolute_10_20_30.pdf')
    plt.rcParams.update(plt.rcParamsDefault)
    plt.close('all')

def plot_2aspect_scatter_times_vs_number_of_subnets_all_algos(net_kw_subnet_generator_list,savename,min_number=0):
    nl_mesu_scatter = [[],[]]
    a_mesu_scatter = [[],[]]
    aggregated_scatter = [[],[]]
    for net_kw_subnet_generator in net_kw_subnet_generator_list:
        res_dict = load_benchmark_result_by_subnet_size_all_algos(net_kw_subnet_generator)
        for nnodes in res_dict:
            for nlayers_in_second_aspect in res_dict[nnodes]:
                for nlayers_in_first_aspect in res_dict[nnodes][nlayers_in_second_aspect]:
                    for subnet_size in res_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect]:
                        alg_res = res_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect][subnet_size]
                        if alg_res['nl-mesu_number'] > min_number or alg_res['a-mesu_number'] > min_number or alg_res['aggregated_number'] > min_number:
                            nl_mesu_scatter[0].append(alg_res['nl-mesu_number'])
                            nl_mesu_scatter[1].append(alg_res['nl-mesu_time'])
                            a_mesu_scatter[0].append(alg_res['a-mesu_number'])
                            a_mesu_scatter[1].append(alg_res['a-mesu_time'])
                            aggregated_scatter[0].append(alg_res['nl-mesu_number'])
                            aggregated_scatter[1].append(alg_res['aggregated_time'])
                        #if alg_res['nl-mesu_number'] == 1:
                        #    print(nnodes,nlayers_in_first_aspect,nlayers_in_second_aspect,subnet_size)
    fig,ax = plt.subplots(figsize=(0.7*6.4,0.7*4.8))
    # fitlines
    fit_x_1 = np.logspace(2,4,num=50,base=10)
    fit_x_2 = np.logspace(3,5,num=50,base=10)
    fit_y_1 = []
    fit_y_2 =  []
    for x in fit_x_1:
        fit_y_1.append(1 * (x**1.2))
    for x in fit_x_2:
        fit_y_2.append(10**-5 * (x**1))
    ax.plot(fit_x_1,fit_y_1,color='k',linestyle='-.',zorder=-1)
    plt.annotate('',xy=(200,700),xytext=(140,9000),arrowprops=dict(arrowstyle='simple,tail_width=0.05,head_width=0.5',facecolor='black'))
    plt.text(100,10000,r'$t = n_{subnets}^{1.2}$',fontsize=10)
    #plt.text(800,20000,r'$t = 10^{-1.5} \times n_{subnets}^{1.2}$',fontsize=10)
    ax.plot(fit_x_2,fit_y_2,color='k',linestyle='-.',zorder=-1)
    plt.annotate('',xy=(8000,0.06),xytext=(6000,0.011),arrowprops=dict(arrowstyle='simple,tail_width=0.05,head_width=0.5',facecolor='black'))
    plt.text(1500,0.005,r'$t = 10^{-5} \times n_{subnets}$',fontsize=10)
    # actual data
    ax.scatter(aggregated_scatter[0],aggregated_scatter[1],facecolors='none',edgecolors='#7570b3',marker='o',s=20,label='ad-esu',linewidth=1)
    ax.scatter(nl_mesu_scatter[0],nl_mesu_scatter[1],color='#d95f02',marker='x',label='nlse',s=20,linewidth=1)
    ax.scatter(a_mesu_scatter[0],a_mesu_scatter[1],color='#1b9e77',marker='+',label='elsse',s=40,linewidth=1)
    #print(a_mesu_scatter[0] == nl_mesu_scatter[0])
    #print(a_mesu_scatter[1] == nl_mesu_scatter[1])
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Number of subnetworks',size=12)
    ax.set_ylabel(r'$t$ (s)',size=12)
    ax.legend(loc=(0.74,0.005),fontsize=10)
    fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
    if savename:
        fig.savefig(savename)
    else:
        return fig,ax

def plot_geo_er_scatter_convenience_script():
    nnodes = [100,200,300]
    #nnodes = [(100)]
    mean_degs = [(1.5,1.5),(2,1),(1,2)]
    #mean_degs = [(1.5,1.5)]
    min_number = 2 # smallest number of subnetworks to be included
    net_kw_subnet_generator_list = []
    for nn in nnodes:
        for d in mean_degs:
            net_kw_subnet_generator_list.append(make_geo_er_mixed_mlayer_generator_param_mean_degrees(nn,d[0],d[0],d[1],d[1]))
    #net_kw_subnet_generator_list = [make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,1.5,1.5,1.5,1.5),make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,2,2,1,1),make_geo_er_mixed_mlayer_generator_param_mean_degrees(100,1,1,2,2)]
    savename = 'cpp_benchmark_figures/geo_er_mixed_param_mean_deg_all_scatter.pdf'
    plot_2aspect_scatter_times_vs_number_of_subnets_all_algos(net_kw_subnet_generator_list,savename,min_number)
    plt.close('all')

def find_subnetwork_numbers_in_range_2aspect(net_kw_subnet_generator,lower_bound=-1,upper_bound=2):
    # subnetwork numbers between (lower_bound,upper_bound) noninclusive
    # uses nl-mesu result
    res_dict,outputfile_dict = load_benchmark_result_by_subnet_size_all_algos(net_kw_subnet_generator,return_outputfiles=True)
    for nnodes in res_dict:
        for nlayers_in_second_aspect in res_dict[nnodes]:
            for nlayers_in_first_aspect in res_dict[nnodes][nlayers_in_second_aspect]:
                for subnet_size in res_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect]:
                    alg_res = res_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect][subnet_size]
                    if alg_res['nl-mesu_number'] < upper_bound and alg_res['nl-mesu_number'] > lower_bound:
                        files = outputfile_dict[nnodes][nlayers_in_second_aspect][nlayers_in_first_aspect][subnet_size]
                        yield (nnodes,nlayers_in_first_aspect,nlayers_in_second_aspect,subnet_size,alg_res['nl-mesu_number'],alg_res['nl-mesu_time'],files)

def find_anomalous_subnetwork_numbers_geo_er_mixed(lower_bound=-1,upper_bound=2,return_files=True):
    nnodes = [100,200,300]
    mean_degs = [(1.5,1.5),(2,1),(1,2)]
    files = []
    for nn in nnodes:
        for d in mean_degs:
            net_kw_subnet_generator = make_geo_er_mixed_mlayer_generator_param_mean_degrees(nn,d[0],d[0],d[1],d[1])
            print(nn,d)
            for anomaly in find_subnetwork_numbers_in_range_2aspect(net_kw_subnet_generator,lower_bound,upper_bound):
                print(anomaly[:-1])
                files.append(anomaly[-1])
    if return_files:
        return files
