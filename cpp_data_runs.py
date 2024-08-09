#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import helpers
import os
import pymnet as pn
import collections
import itertools
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import random

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

def plot_group_of_lines_xlabels(shared_x_axis,x_axis_labels,individual_y_axes,formats,line_labels,x_axis_label,y_axis_label,title,plot_to_ax=None,colors=None,linewidths=None):
    if plot_to_ax is None:
        fig,ax = plt.subplots()
    else:
        ax = plot_to_ax
        fig = ax.get_figure()
    for ii in range(len(individual_y_axes)):
        if colors and linewidths:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii],color=colors[ii],linewidth=linewidths[ii],alpha=1)
        elif colors:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii],color=colors[ii],alpha=1)
        elif linewidths:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii],linewidth=linewidths[ii],alpha=1)
        else:
            ax.plot(shared_x_axis,individual_y_axes[ii],formats[ii],label=line_labels[ii])
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(y_axis_label)
    ax.set_title(title)
    ax.set_xticks(ticks=shared_x_axis,labels=x_axis_labels)
    return fig,ax

def load_all_algos_results(base_filename):
    # for fb, airline, twitter style result saving (all in one folder)
    outputfiles = []
    # nl-mesu and a-mesu
    outputfiles.append(base_filename+'_both')
    outputfiles.append(base_filename+'_aggregated')
    # defaults to None
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
    return res_dict

def convert_ppi_nets():
    for name in ids:
        if name == 'sacchcere':
            input_filename = 'multiplex_pp_data/SacchCere_Multiplex_Genetic/Dataset/'+name+'_genetic_multiplex.edges'
        elif name == 'sacchpomb':
            input_filename = 'multiplex_pp_data/SacchPomb_Multiplex_Genetic/Dataset/'+name+'_genetic_multiplex.edges'
        else:
            input_filename = 'multiplex_pp_data/'+name.capitalize()+'_Multiplex_Genetic/Dataset/'+name+'_genetic_multiplex.edges'
        output_filename = 'reformatted_ppi_data/'+name+'.edges'
        helpers.reformat_ppi_file(input_filename,output_filename)
    

def run_one_ppi_net_aggregated(name):
    net_folder = 'reformatted_ppi_data/'
    result_folder = 'cpp_results_aggregated_algo/'
    inputfile = net_folder+name+'.edges'
    n_aspects = 1
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    for size in sizes:
        outputfile = result_folder + name + "_(" + str(size[0]) + "," + str(size[1]) + ").txt"
        call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(size).replace(' ','').strip('()') + "'"
        call_str = call_str + ' time '
        call_str = call_str + ' aggregated'
        if not os.path.exists(outputfile):
            os.system(call_str)

def get_ppi_times():
    both_result_folder = 'cpp_results/'
    agg_result_folder = 'cpp_results_aggregated_algo/'
    full_res_dict = {}
    for org_id in ids:
        full_res_dict[org_id] = dict()
        for size in sizes:
            outputfiles = []
            for result_folder in [both_result_folder,agg_result_folder]:
                outputfiles.append(result_folder + org_id + "_(" + str(size[0]) + "," + str(size[1]) + ").txt")
            # defaults to None
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
            full_res_dict[org_id][size] = res_dict
    return full_res_dict

##### FB network

def make_fb_network(timescale='days'):
    # fb messaging network
    # does not include self-loops
    save_folder = 'cpp_fb_networks/'
    net=pn.MultilayerNetwork(fullyInterconnected=False,aspects=1)
    if timescale == 'hours':
        seconds_divisor = 60*60
    elif timescale == 'days':
        seconds_divisor = 60*60*24
    elif timescale == 'weeks':
        seconds_divisor = 60*60*24*7
    elif timescale == 'threedays':
        seconds_divisor = 60*60*24*3
    # intralayer links
    with open("fb_data/fb-messages.edges",'r') as f:
        for line in f:
            fr,to,time=line.split(',')
            fr,to = int(fr),int(to)
            time= int(float(time))
            time_layer=int((time-1080090715)/seconds_divisor)
            if fr != to:
                net[fr,to,time_layer,time_layer]=1
    # interlayer links
    # forward to next occurrence of node
    node_to_layers=collections.defaultdict(list)
    for node,layer in net.iter_node_layers():
        node_to_layers[node].append(layer)
    for node,layers in node_to_layers.items():
        layers=sorted(layers)
        for i in range(len(layers)):
            if i!=0:
                net[node,node,layers[i-1],layers[i]]=1
    # save in cpp format
    helpers.save_edgelist_cpp_format(net,savename=save_folder+timescale)

def run_fb_network(timescale='days',algo='both'):
    result_folder = 'cpp_fb_results/'
    n_aspects = 1
    inputfile = 'cpp_fb_networks/'+timescale
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    for size in sizes:
        outputfile = result_folder+timescale+'_(' + str(size[0]) + ',' + str(size[1]) + ')_' + algo
        call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(size).replace(' ','').strip('()') + "'"
        call_str = call_str + ' time '
        call_str = call_str + ' ' + algo
        if not os.path.exists(outputfile):
            os.system(call_str)

def get_fb_net_times():
    algos = ['both','aggregated']
    timescales = ['hours','days','threedays','weeks']
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    result_folder = 'cpp_fb_results/'
    return_dict = dict()
    for ts in timescales:
        return_dict[ts] = dict()
        for size in sizes:
            return_dict[ts][size] = dict()
            nl_found = False
            a_found = False
            agg_found = False
            for algo in algos:
                filename = result_folder + ts + '_(' + str(size[0]) + ',' + str(size[1]) + ')_' + algo
                try:
                    with open(filename,'r') as f:
                        file_is_empty = True
                        for line in f:
                            data = line.split()
                            if data[0] == 'nl-mesu':
                                nl_mesu_time = float(data[1])
                                nl_mesu_number = int(data[2])
                                return_dict[ts][size]['nl-mesu'] = (nl_mesu_time,nl_mesu_number)
                                nl_found = True
                            elif data[0] == 'a-mesu':
                                a_mesu_time = float(data[1])
                                a_mesu_number = int(data[2])
                                return_dict[ts][size]['a-mesu'] = (a_mesu_time,a_mesu_number)
                                a_found = True
                            elif data[0] == 'aggregated':
                                aggregated_esu_time = float(data[1])
                                aggregated_esu_number = int(data[2])
                                return_dict[ts][size]['aggregated'] = (aggregated_esu_time,aggregated_esu_number)
                                agg_found = True
                            file_is_empty = False
                        if file_is_empty:
                            raise Exception # go to except block
                        #assert nl_mesu_number == a_mesu_number and nl_mesu_number == aggregated_esu_number
                except:
                    pass
            if not nl_found:
                return_dict[ts][size]['nl-mesu'] = (None,None)
            if not a_found:
                return_dict[ts][size]['a-mesu'] = (None,None)
            if not agg_found:
                return_dict[ts][size]['aggregated'] = (None,None)
    return return_dict

def get_fb_times():
    result_folder = 'cpp_fb_results/'
    res_dict = {}
    for agg_level in ['hours','days','threedays','weeks']:
        savename = agg_level
        res_dict[agg_level] = dict()
        for size in sizes:
            base_filename = result_folder+savename+'_(' + str(size[0]) + ',' + str(size[1]) + ')'
            res_dict[agg_level][size] = load_all_algos_results(base_filename)
    return res_dict

def plot_fb_net_times():
    timescales = ['hours','days','threedays','weeks']
    #sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    sizes = ((2,2),(2,3),(3,2),(3,3),(4,2),(4,3))
    algos = ['nl-mesu','a-mesu','aggregated']
    fb_net_times = get_fb_net_times()
    colors = {'hours' : "#b30000", 'days' : "#e34a33", 'threedays' : "#fc8d59", 'weeks' : "#fdbb84"}
    formats = {'nl-mesu' : '-', 'a-mesu' : '--', 'aggregated' : ':'}
    #colors = ["#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"]
    all_y_vals = []
    all_formats = []
    all_colors = []
    all_labels = []
    for ts in timescales:
        for algo in algos:
            y_vals = []
            for size in sizes:
                y_vals.append(fb_net_times[ts][size][algo][0])
            all_y_vals.append(y_vals)
            all_formats.append(formats[algo])
            all_colors.append(colors[ts])
            all_labels.append((ts,size,algo))
    fig,ax = plt.subplots(1,1,figsize=(0.79*6.4,0.79*4.8))
    shared_x_axis = list(range(len(sizes)))
    xlabs = [str(size) for size in sizes]
    plot_group_of_lines_xlabels(shared_x_axis = shared_x_axis, x_axis_labels = xlabs, individual_y_axes = all_y_vals, formats = all_formats, line_labels = all_labels, x_axis_label = 'Subnetwork size', y_axis_label = 'Running time (s)', title = '', plot_to_ax = ax, colors = all_colors)
    ax.set_yscale('log')
    # legends
    
    lines2 = [Line2D([0], [0], color='k', linewidth=1, linestyle=ls) for ls in ['-','--','-.']]
    labels2 = ['nlse','elsse','ad-esu']
    legend2 = plt.legend(lines2, labels2, loc=8)
    
    colors = ['#b30000', '#e34a33', '#fc8d59', '#fdbb84']
    lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
    labels = ['Hour', 'Day', '3 days', 'Week']
    plt.legend(lines, labels, loc=4)
    
    ax.add_artist(legend2)
    fig.savefig('cpp_fb_figures/fb_absolute_running_times.pdf')
            
##### Airline network

def make_airline_network(n_layers='all'):
    save_folder = 'cpp_airline_networks/'
    savename = '_airlines'
    net=pn.MultilayerNetwork(fullyInterconnected=False,aspects=1)
    curr_layer = -1 # 0 will be first layer
    with open('airline_data/network.txt','r') as f:
        for line in f:
            line_split = line.split()
            if len(line_split) == 1:
                curr_layer = curr_layer + 1
            elif len(line_split) == 0: # empty line between layer header and edges
                pass
            else: # edge information
                curr_node = int(line_split[0])
                deg = int(line_split[1])
                for neighbor in line_split[2:]:
                    net[curr_node,int(neighbor),curr_layer,curr_layer] = 1
                    #print('intralayer edge')
            node_to_layers=collections.defaultdict(list)
    # interlayer links
    # all connected to all
    ils = 0
    node_to_layers=collections.defaultdict(list)
    for node,layer in net.iter_node_layers():
        node_to_layers[node].append(layer)
    for node,layers in node_to_layers.items():
        layers=sorted(layers)
        for layer_pair in itertools.combinations(layers,2):
            net[node,node,layer_pair[0],layer_pair[1]] = 1
            #print('interlayer edge')
            ils = ils+1
    print(ils)
    print(len(net.edges))
    # save in cpp format
    print(len(list(net.iter_layers())))
    #print(node_to_layers)
    if n_layers == 'all':
        savename = 'all' + savename
        helpers.save_edgelist_cpp_format(net,savename=save_folder+savename)
    else:
        top_n_net = get_top_layers_with_most_edges(net,n_layers)
        savename = str(n_layers) + savename
        print(len(list(top_n_net.iter_layers())))
        print(len(top_n_net.edges))
        helpers.save_edgelist_cpp_format(top_n_net,savename=save_folder+savename)

def get_top_layers_with_most_edges(net,n=10):
    # connects nodes on all layers to their counterparts on other layers
    new_net=pn.MultilayerNetwork(fullyInterconnected=False,aspects=1)
    edges_per_layer = dict()
    for ll in net.iter_layers():
        edges_per_layer[ll] = 0
    for ee in net.edges:
        if ee[2] == ee[3]: # intralayer edge
            edges_per_layer[ee[2]] = edges_per_layer[ee[2]] + 1
    for layer,n_edges in sorted(edges_per_layer.items(),key=lambda x:x[1],reverse=True)[0:n]:
        for node in net.iter_nodes(layer=layer):
            for neigh in net[(node,layer)]:
                if layer == neigh[1]: # neighbor on the same layer
                    new_net[node,neigh[0],layer,layer] = 1
    # add interlayer edges
    node_to_layers=collections.defaultdict(list)
    for node,layer in new_net.iter_node_layers():
        node_to_layers[node].append(layer)
    for node,layers in node_to_layers.items():
        layers=sorted(layers)
        for layer_pair in itertools.combinations(layers,2):
            new_net[node,node,layer_pair[0],layer_pair[1]] = 1
    return new_net

def run_airline_network(savename='all_airlines',algo='both'):
    result_folder = 'cpp_airline_results/'
    n_aspects = 1
    inputfile = 'cpp_airline_networks/'+savename
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    for size in sizes:
        outputfile = result_folder+savename+'_(' + str(size[0]) + ',' + str(size[1]) + ')_' + algo
        call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(size).replace(' ','').strip('()') + "'"
        call_str = call_str + ' time '
        call_str = call_str + ' ' + algo
        if not os.path.exists(outputfile):
            os.system(call_str)

def get_airline_times():
    result_folder = 'cpp_airline_results/'
    res_dict = {}
    for nairlines in ['10','20','all']:
        savename = nairlines+'_airlines'
        res_dict[nairlines] = dict()
        for size in sizes:
            base_filename = result_folder+savename+'_(' + str(size[0]) + ',' + str(size[1]) + ')'
            res_dict[nairlines][size] = load_all_algos_results(base_filename)
    return res_dict

##### twitter network

def make_file_ids():
    single_hashtags = set()
    topics = set()
    for fname in sorted(os.listdir('twitter_data/')):
        if fname[-9:] == '.edgelist':
            # remove p1,p2,p3
            fname = fname.split('_')[0]
            # check for first letter being uppercase
            if fname[0].isupper():
                topics.add(fname)
            else:
                single_hashtags.add(fname)
    return {'single_hashtags' : sorted(list(single_hashtags)), 'topics' : sorted(list(topics))}

def make_twitter_network(level='single_hashtags',linkage='all'):
    net = pn.MultilayerNetwork(fullyInterconnected=False,aspects=1)
    data_dir = 'twitter_data/'
    save_folder = 'cpp_twitter_networks/'
    layer = -1 # layers as integers
    file_id_to_layers = {}
    for file_id in make_file_ids()[level]:
        layer = layer + 1
        file_id_to_layers[file_id] = layer
        for time in ['_p1','_p2','_p3']:
            try:
                with open(data_dir+file_id+time+'.edgelist','r') as f:
                    for line in f:
                        if not line[0:3] == 'ret': # eliminate header line
                            node1,node2,weight = line.split(',')
                            node1 = int(node1)
                            node2 = int(node2)
                            net[node1,node2,layer,layer] = 1
            except:
                    pass
    # make node_to_layers
    node_to_layers=collections.defaultdict(list)
    for node,layer in net.iter_node_layers():
        node_to_layers[node].append(layer)
    # make layer_to_nodes
    layer_to_nodes = collections.defaultdict(set)
    for ll in net.iter_layers():
        layer_to_nodes[ll] = set(net.iter_nodes(layer=ll))
    # make interlayer edges according to linkage
    if linkage == 'all':
        # interlayer edges, all to all
        for node,layers in node_to_layers.items():
            layers=sorted(layers)
            for layer_pair in itertools.combinations(layers,2):
                net[node,node,layer_pair[0],layer_pair[1]] = 1
    elif linkage == 'max_jaccard_linkage':
        for ll1 in net.iter_layers():
            # if multiple max, choose a random one
            max_jaccards = [0]
            max_jaccard_layers = []
            for ll2 in net.iter_layers():
                if ll1 != ll2:
                    curr_jaccard = len(layer_to_nodes[ll1].intersection(layer_to_nodes[ll2]))/len(layer_to_nodes[ll1].union(layer_to_nodes[ll2]))
                    if curr_jaccard > max_jaccards[0]:
                        max_jaccards = [curr_jaccard]
                        max_jaccard_layers = [ll2]
                    elif curr_jaccard == max_jaccards[0]:
                        max_jaccards.append(curr_jaccard)
                        max_jaccard_layers.append(ll2)
            # forward to one, max jaccard
            target_layer = random.choice(max_jaccard_layers)
            for shared_node in layer_to_nodes[ll1].intersection(layer_to_nodes[target_layer]):
                net[shared_node,shared_node,ll1,target_layer] = 1
    elif linkage == 'max_NMI_linkage':
        # get max NMI values from Chen, T. H. Y., Salloum, A., Gronow, A., Ylä-Anttila, T., & Kivelä, M. (2021). Polarization of climate politics results from partisan sorting: Evidence from Finnish Twittersphere. Global Environmental Change, 71, 102348.
        # using NMI from pre-election period
        # and set forward linkage source : target(s)
        max_nmis = {'CLIMATE' : ['LEFT'],
         'IMMIGRATION' : ['FINNS'],
         'SOCIALSECURITY' : ['CENTRE'],
         'ECONOMICPOLICY' : ['SDP','FINNS'], # two equal ones
         'EDUCATION' : ['SDP'],
         'SDP' : ['LEFT'],
         'FINNS' : ['IMMIGRATION'],
         'NATIONAL' : ['SOCIALSECURITY'],
         'CENTRE' : ['SOCIALSECURITY'],
         'GREEN' : ['CLIMATE','LEFT'], # two equal ones
         'LEFT' : ['IMMIGRATION']}
        for source_layer in max_nmis:
            for target_layer in max_nmis[source_layer]:
                for shared_node in layer_to_nodes[file_id_to_layers[source_layer]].intersection(layer_to_nodes[file_id_to_layers[target_layer]]):
                    net[shared_node,shared_node,file_id_to_layers[source_layer],file_id_to_layers[target_layer]] = 1
    elif linkage == 'none':
        # no interlayer edges for comparison
        pass
    else:
        raise NotImplementedError()
    savename = level+'_'+linkage
    helpers.save_edgelist_cpp_format(net,save_folder+savename)

def run_twitter_network(level='single_hasthatgs',linkage='all',algo='both'):
    # TODO NB !!!!
    # uses the temp file with buffer flushing endl !!!!!
    result_folder = 'cpp_twitter_results/'
    n_aspects = 1
    identifier = level+'_'+linkage
    inputfile = 'cpp_twitter_networks/'+identifier
    sizes = ((2,2),(3,2),(2,3),(3,3),(4,2),(4,3))
    for size in sizes:
        outputfile = result_folder+identifier+'_(' + str(size[0]) + ',' + str(size[1]) + ')_' + algo
        #call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(size).replace(' ','').strip('()') + "'"
        call_str = './mesu_' + str(n_aspects) + '_TEMP_FOR_ENDL.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(size).replace(' ','').strip('()') + "'"
        call_str = call_str + ' time '
        call_str = call_str + ' ' + algo
        if not os.path.exists(outputfile):
            os.system(call_str)

def get_twitter_times():
    result_folder = 'cpp_twitter_results/'
    res_dict = {}
    for net_type in ['single_hashtags','topics','single_hashtags_max_jaccard_linkage','topics_max_NMI_linkage']:
        savename = net_type
        res_dict[net_type] = dict()
        for size in sizes:
            base_filename = result_folder+savename+'_(' + str(size[0]) + ',' + str(size[1]) + ')'
            res_dict[net_type][size] = load_all_algos_results(base_filename)
    return res_dict

##### full results aggregation

def get_all_results():
    all_res = {}
    all_res['ppi'] = get_ppi_times()
    all_res['fb'] = get_fb_times()
    all_res['airline'] = get_airline_times()
    all_res['twitter'] = get_twitter_times()
    return all_res

def filter_all_results(all_res,filter_type='two_of_each'):
    if filter_type == 'two_of_each':
        filtered_allowed = {'ppi' : ['celegans','arabidopsis'],
                            'fb' : ['hours','weeks'],
                            'airline' : ['20','all'],
                            'twitter' : ['single_hashtags_max_jaccard_linkage','topics_max_NMI_linkage']}
    elif filter_type == 'ppi_air':
        filtered_allowed = {'ppi' : ['celegans','arabidopsis'],
                            'fb' : [],
                            'airline' : ['20','all'],
                            'twitter' : []}
    elif filter_type == 'fb_twitter':
        filtered_allowed = {'ppi' : [],
                            'fb' : ['hours','weeks'],
                            'airline' : [],
                            'twitter' : ['single_hashtags_max_jaccard_linkage','topics_max_NMI_linkage']}
    filtered_res = {}
    for dataset in all_res:
        filtered_res[dataset] = dict()
        for net_type in all_res[dataset]:
            if net_type in filtered_allowed[dataset]:
                filtered_res[dataset][net_type] = all_res[dataset][net_type]
    return filtered_res

def convert_results_to_lists(all_res,sizes=[(2,2),(2,3),(3,2),(3,3),(4,2),(4,3)]):
    list_dict= dict()
    for dataset in all_res:
        list_dict[dataset] = dict()
        for net_type in all_res[dataset]:
            list_dict[dataset][net_type] = dict()
            for algo in ['nl-mesu','a-mesu','aggregated']:
                list_dict[dataset][net_type][algo] = []
                for size in sizes:
                        list_dict[dataset][net_type][algo].append(all_res[dataset][net_type][size][algo+'_time'])
    return list_dict

def format_results_to_plotting(all_res_list_format,data_colors):
    formats_dict = {'nl-mesu' : '-', 'a-mesu' : '--', 'aggregated' : ':'}
    line_widths_dict = {'nl-mesu' : 2, 'a-mesu' : 2.5, 'aggregated' : 3}
    individual_y_axes = []
    formats = []
    line_labels = []
    colors = []
    linewidths = []
    for dataset in all_res_list_format:
        for net_type in all_res_list_format[dataset]:
            for algo in ['nl-mesu','a-mesu','aggregated']:
                individual_y_axes.append(all_res_list_format[dataset][net_type][algo])
                formats.append(formats_dict[algo])
                line_labels.append(dataset+'_'+net_type+'_'+algo)
                colors.append(data_colors[dataset][net_type])
                linewidths.append(line_widths_dict[algo])
    return individual_y_axes,formats,line_labels,colors,linewidths

def plot_all_results(included='all',plot_to_ax=None):
    plt.rcParams.update({'font.size':12, 'legend.fontsize': 12,'legend.handlelength': 3.5,'legend.loc':'upper center','legend.columnspacing': 0.7,'legend.handletextpad': 0.3,'lines.linewidth':3})
    # final publication image:
    # ppi: celegans, arabidopsis
    # fb: hours (comparison mostly nl-mesu and a-mesu), weeks
    # airline: top 20 , all
    # twitter: hashtags Jaccard, topics NMI
    all_res = get_all_results()
    super_colors = {'ppi' : '#E69F00', 'fb' : '#000000', 'airline' : '#009E73', 'twitter' : '#CC79A7'}
    individual_colors = {'ppi' : {'celegans':'#FFC033','arabidopsis':'#E69F00'},
                         'fb' : {'hours':'#000000','weeks':'#757575'},
                         'airline' : {'20':'#33FFC7','all':'#00996F'},
                         'twitter' : {'single_hashtags_max_jaccard_linkage':'#94386B','topics_max_NMI_linkage':'#D590B6'}}
    if included == 'all':
        individual_colors = collections.defaultdict(dict)
        for dataset in all_res:
            for net_type in all_res[dataset]:
                individual_colors[dataset][net_type] = super_colors[dataset]
    sizes = [(2,2),(2,3),(3,2),(3,3),(4,2),(4,3)]
    if included == 'final_publication':
        all_res = filter_all_results(all_res)
    if included == 'ppi_air':
        all_res = filter_all_results(all_res,included)
    if included == 'fb_twitter':
        all_res = filter_all_results(all_res,included)
    all_res_list_format = convert_results_to_lists(all_res,sizes)
    individual_y_axes,formats,line_labels,colors,linewidths = format_results_to_plotting(all_res_list_format, individual_colors)
    shared_x_axis = list(range(len(sizes)))
    x_axis_labels = [''.join(str(x).split()) for x in sizes]
    if not plot_to_ax:
        fig,ax = plt.subplots(1,1)
    else:
        ax = plot_to_ax
    if plot_to_ax:
        plot_group_of_lines_xlabels(shared_x_axis, x_axis_labels, individual_y_axes, formats, line_labels, x_axis_label=None, y_axis_label=None, title='', plot_to_ax=ax, colors=colors, linewidths=linewidths)
        ax.set_yscale('log')
        if included == 'fb_twitter':
            ax.scatter([0],all_res['fb']['hours'][(2,2)]['aggregated_time'],marker=',',facecolors=individual_colors['fb']['hours'],edgecolors='none',s=10)
            ax.text(-0.18,1.5*10**5,'FB hourly AD-ESU')
            plt.annotate('',xy=(0.06,all_res['fb']['hours'][(2,2)]['aggregated_time']+3000),xytext=(0.6,1.2*10**5),arrowprops=dict(arrowstyle='simple,tail_width=0.03,head_width=0.4',facecolor='black'))
    if not plot_to_ax:
        plot_group_of_lines_xlabels(shared_x_axis, x_axis_labels, individual_y_axes, formats, line_labels, x_axis_label='Subnetwork size', y_axis_label='$t$ (s)', title='', plot_to_ax=ax, colors=colors)
        ax.set_yscale('log')
        fig.subplots_adjust(top=0.934,right=0.99,left=0.115,bottom=0.095)
        if included == 'all':
            fig.text(0.13,0.95,'Dataset:',fontsize=10)
            addition = 0.12
            for dataset in ['ppi','fb','airline','twitter']:
                fig.text(0.13+addition,0.95,dataset,color=super_colors[dataset],fontsize=10)
                addition = addition + 0.02 + 0.01*len(dataset)
        elif included == 'final_publication':
            fig.text(0.001,0.944,'C. elegans',color=individual_colors['ppi']['celegans'],fontsize=10,weight='bold',style='italic')
            fig.text(0.14,0.944,'A. thaliana',color=individual_colors['ppi']['arabidopsis'],fontsize=10,weight='bold',style='italic')
            fig.text(0.285,0.944,'FB hourly',color=individual_colors['fb']['hours'],fontsize=10,weight='bold')
            fig.text(0.41,0.944,'FB weekly',color=individual_colors['fb']['weeks'],fontsize=10,weight='bold')
            fig.text(0.545,0.944,'Air all',color=individual_colors['airline']['all'],fontsize=10,weight='bold')
            fig.text(0.627,0.944,'Air 20',color=individual_colors['airline']['20'],fontsize=10,weight='bold')
            fig.text(0.71,0.944,'Twitter #',color=individual_colors['twitter']['single_hashtags_max_jaccard_linkage'],fontsize=10,weight='bold')
            fig.text(0.83,0.944,'Twitter topics',color=individual_colors['twitter']['topics_max_NMI_linkage'],fontsize=10,weight='bold')
        lines2 = [Line2D([0], [1], color='k', linewidth=2, linestyle=ls) for ls in ['-','--',':']]
        labels2 = ['nlse','elsse','ad-esu']
        legend2 = plt.legend(lines2, labels2, loc=8, ncols=3, bbox_to_anchor=(0.5,1),frameon=False,fancybox=False,shadow=False)
        # add individual point for aggregated FB hours
        ax.scatter([0],all_res['fb']['hours'][(2,2)]['aggregated_time'],marker=',',facecolors=individual_colors['fb']['hours'],edgecolors='none',s=10)
        fig.text(0.13,0.79,'FB hourly AD-ESU')
        plt.annotate('',xy=(0.06,all_res['fb']['hours'][(2,2)]['aggregated_time']+3000),xytext=(0.5,80000),arrowprops=dict(arrowstyle='simple,tail_width=0.03,head_width=0.4',facecolor='black'))
        plt.savefig('cpp_figures/absolute_running_times_extradata_'+included+'.pdf')
        plt.close('all')
    plt.rcParams.update(plt.rcParamsDefault)

def plot_all_results_split():
    plt.rcParams.update({'font.size':12, 'legend.fontsize': 12,'legend.handlelength': 3.5,'legend.loc':'upper center','legend.columnspacing': 0.7,'legend.handletextpad': 0.3,'lines.linewidth':3})
    individual_colors = {'ppi' : {'celegans':'#FFC033','arabidopsis':'#E69F00'},
                         'fb' : {'hours':'#000000','weeks':'#757575'},
                         'airline' : {'20':'#33FFC7','all':'#00996F'},
                         'twitter' : {'single_hashtags_max_jaccard_linkage':'#94386B','topics_max_NMI_linkage':'#D590B6'}}
    # multiply width with 1.1 for width 0.9 -> 0.99
    fig,ax = plt.subplots(1,2,sharey='all',figsize=(1.1*1.8*0.79*6.4,0.79*4.8))
    plot_all_results('ppi_air',ax[0])
    plot_all_results('fb_twitter',ax[1])
    fig.subplots_adjust(wspace=0, hspace=0)
    #fig.subplots_adjust(top=0.935,right=0.995,left=0.06,bottom=0.1)
    fig.subplots_adjust(top=0.917,right=0.995,left=0.065,bottom=0.11)
    fig.supxlabel('          Subnetwork size',y=0)
    fig.supylabel(r'           $t$ (s)',x=0)
    # legends
    plt.rcParams.update({'font.size':12, 'legend.fontsize': 12,'legend.handlelength': 3.5,'legend.loc':'upper center','legend.columnspacing': 0.7,'legend.handletextpad': 0.3,'lines.linewidth':3})
    legend_linewidths_by_style = {'-' : 2, '--' : 2.5, ':' : 3}
    lines2 = [Line2D([0], [1], color='k', linewidth=legend_linewidths_by_style[ls], linestyle=ls) for ls in ['-','--',':']]
    # turn on latex and setup for \textsc
    #plt.rc('text', usetex=True)
    #plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage[T1]{fontenc} \usepackage{mathptmx}')
    #plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    #labels2 = [r'$\textbf{NLSE}$',r'$\textbf{\textsc{elsse}}$',r'$\textsc{ad-esu}$']
    labels2 =  ['NLSE','ELSSE','AD-ESU'] # just use all caps and times new roman font
    legend2 = plt.legend(lines2, labels2, loc=8, ncols=3, bbox_to_anchor=(0,0.999),frameon=False,fancybox=False,shadow=False, prop={'family':'Times New Roman'})
    # turn off latex for the rest of the text
    #plt.rc('text', usetex=False)
    y_addition = -0.018
    fig.text(0.12,y_addition+0.944,'C. elegans',color=individual_colors['ppi']['celegans'],fontsize=10,weight='bold',style='italic')
    fig.text(0.22,y_addition+0.944,'A. thaliana',color=individual_colors['ppi']['arabidopsis'],fontsize=10,weight='bold',style='italic')
    fig.text(0.57,y_addition+0.944,'FB hourly',color=individual_colors['fb']['hours'],fontsize=10,weight='bold')
    fig.text(0.66,y_addition+0.944,'FB weekly',color=individual_colors['fb']['weeks'],fontsize=10,weight='bold')
    fig.text(0.326,y_addition+0.944,'Air all',color=individual_colors['airline']['all'],fontsize=10,weight='bold')
    fig.text(0.39,y_addition+0.944,'Air 20',color=individual_colors['airline']['20'],fontsize=10,weight='bold')
    fig.text(0.755,y_addition+0.944,'Twitter #',color=individual_colors['twitter']['single_hashtags_max_jaccard_linkage'],fontsize=10,weight='bold')
    fig.text(0.84,y_addition+0.944,'Twitter topics',color=individual_colors['twitter']['topics_max_NMI_linkage'],fontsize=10,weight='bold')
    fig.savefig('cpp_figures/absolute_running_times_extradata_final_publication_split_v2.pdf')
    plt.close('all')
    plt.rcParams.update(plt.rcParamsDefault)


def format_results_to_scatter(all_res):
    count = 0
    nl_mesu_scatter = [[],[]]
    a_mesu_scatter = [[],[]]
    aggregated_scatter = [[],[]]
    for dataset in all_res:
        for net in all_res[dataset]:
            for subnet_size in all_res[dataset][net]:
                count = count + 1
                nl_mesu_scatter[0].append(all_res[dataset][net][subnet_size]['nl-mesu_number'])
                nl_mesu_scatter[1].append(all_res[dataset][net][subnet_size]['nl-mesu_time'])
                a_mesu_scatter[0].append(all_res[dataset][net][subnet_size]['a-mesu_number'])
                a_mesu_scatter[1].append(all_res[dataset][net][subnet_size]['a-mesu_time'])
                aggregated_scatter[0].append(all_res[dataset][net][subnet_size]['aggregated_number'])
                aggregated_scatter[1].append(all_res[dataset][net][subnet_size]['aggregated_time'])
    #print(count)
    return nl_mesu_scatter,a_mesu_scatter,aggregated_scatter

def plot_all_results_scatter():
    savename = 'cpp_figures/all_data_scatter.pdf'
    all_res = get_all_results()
    nl_mesu_scatter,a_mesu_scatter,aggregated_scatter = format_results_to_scatter(all_res)
    fig,ax = plt.subplots(figsize=(0.7*6.4,0.7*4.8))
    ax.scatter(nl_mesu_scatter[0],nl_mesu_scatter[1],color='darkred',marker='x',label='nlse',linewidth=1)
    ax.scatter(a_mesu_scatter[0],a_mesu_scatter[1],color='darkgreen',marker='+',label='elsse',linewidth=1)
    ax.scatter(aggregated_scatter[0],aggregated_scatter[1],facecolors='none',edgecolors='darkblue',marker='o',s=20,label='ad-esu',linewidth=1)
    print(a_mesu_scatter[0] == nl_mesu_scatter[0])
    print(a_mesu_scatter[1] == nl_mesu_scatter[1])
    print(a_mesu_scatter[0] == aggregated_scatter[0])
    print(len(a_mesu_scatter[0]))
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Number of subnetworks',size=12)
    ax.set_ylabel(r'$t$ (s)',size=12)
    ax.legend()
    fig.subplots_adjust(top=0.99,right=0.99,left=0.15,bottom=0.13)
    if savename:
        fig.savefig(savename)
    else:
        return fig,ax

def find_min_max_subnet_amounts():
    maxi = 0
    mini = 1000000000000000000
    all_res = filter_all_results(get_all_results())
    print(all_res.keys())
    for i in all_res:
        for j in all_res[i]:
            for k in all_res[i][j]:
                n1 = all_res[i][j][k]['nl-mesu_number']
                n2 = all_res[i][j][k]['a-mesu_number']
                n3 = all_res[i][j][k]['aggregated_number']
                if n1:
                    n = n1
                elif n2:
                    n = n2
                elif n3:
                    n = n3
                if n and n > maxi:
                    maxi = n
                    maxa = (i,j,k)
                if n and n < mini:
                    mini = n
                    mina = (i,j,k)
    print('Max:')
    print(maxi)
    print(maxa)
    print('Min:')
    print(mini)
    print(mina)




















