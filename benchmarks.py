#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import mesu
import time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
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
            ('multiplex_pp_data/Rattus_Multiplex_Genetic/Dataset/rattus_genetic_multiplex.edges',base_p,'rattus')]
#            ('multiplex_pp_data/SacchCere_Multiplex_Genetic/Dataset/sacchcere_genetic_multiplex.edges',0.0001,'sacchcere'),
#            ('multiplex_pp_data/SacchPomb_Multiplex_Genetic/Dataset/sacchpomb_genetic_multiplex.edges',0.0001,'sacchpomb')]
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














