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
    fig,ax = plot_group_of_lines(shared_x_axis, individual_y_axes_nl_mesu, ['-']*len(individual_y_axes_nl_mesu), ids, 'Subnet size', 'Running time (s)', 'Protein-protein interaction multiplex networks',colors=colors)
    plot_group_of_lines(shared_x_axis, individual_y_axes_a_mesu, ['--']*len(individual_y_axes_a_mesu), ids, 'Subnet size', 'Running time (s)', 'Protein-protein interaction multiplex networks',plot_to_ax=ax,colors=colors)
    legend_elements = []
    for ii,name in enumerate(ids):
        legend_elements.append(plt.Line2D([0],[0],marker='o',color=colors[ii],label=name,markerfacecolor=colors[ii],markersize=5,linewidth=0))
    if legend:
        ax.legend(handles=legend_elements,loc='center left',bbox_to_anchor=(1, 0.5))
    ax.set_yscale('log')
    fig.tight_layout()
    fig.savefig(savename)

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
    text_y_multipliers = [0.9,0.8,1.2,1,1,1,1.1,1,1,1,1]
    fig,ax = plt.subplots(figsize=(1.9,4.8))
    colors = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"]
    for ii,ee in enumerate(n_edges):
        ax.scatter(x_values[ii],[ee],c=colors[ii],marker='o')
        plt.text(0.1,ee*text_y_multipliers[ii],ids[ii],ha='left')
    ax.set_xlim([-0.04,0.4])
    ax.set_yscale('log')
    ax.set_ylabel('Number of edges')
    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig('cpp_figures/edge_numbers.pdf')

def make_table_run_times_for_cpp(filename):
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
                curr_line = preamble + " & " + str(size[0]) + ", " + str(size[1]) + " & " + (str(d[name][size][2]) if d[name][size] else "-") + " & " + (str(d[name][size][0]) if d[name][size] else "-") + " & " + (str(d[name][size][1]) if d[name][size] else "-") + r" \\" + "\n"
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

def parse_output_savename(network_inputfilename, subnet_size):
    base_name = network_inputfilename.split('/')[1]
    save_folder = 'cpp_benchmark_results/'
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
    nlayers = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
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
    nlayers = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
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

def plot_mplex_relative_vs_net_size(net_kw_subnet_generator=make_geo_mplex_generator(),savename='cpp_benchmark_figures/geo_mplex.pdf',title='GEO mplex'):
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
    fig,ax = plot_group_of_lines(shared_x_axis,individual_y_axes,formats,line_labels,'Number of layers',r'$t_{a-mesu}/t_{nl-mesu}$',title)
    ax.set_yscale('log')
    ax.legend()
    ax.plot(shared_x_axis,[1]*len(shared_x_axis),'--k')
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
    for kk in sorted(res_dict_with_separate_a_and_nl['nl-mesu'].keys()):
        individual_y_axes.append(res_dict_with_separate_a_and_nl['nl-mesu'][kk])
        formats.append('-')
        line_labels.append(str(kk))
    for kk in sorted(res_dict_with_separate_a_and_nl['a-mesu'].keys()):
        individual_y_axes.append(res_dict_with_separate_a_and_nl['a-mesu'][kk])
        formats.append('--')
        line_labels.append(str(kk))
    fig,ax = plot_group_of_lines(shared_x_axis,individual_y_axes,formats,line_labels,'Number of layers',r'$t$',title)
    ax.set_yscale('log')
    ax.legend()
    #ax.plot(shared_x_axis,[1]*len(shared_x_axis),'--k')
    fig.savefig(savename)

