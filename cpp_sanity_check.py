from helpers import save_edgelist_cpp_format
import os
import sys
import itertools
import pymnet
import pickle

def make_all_small_nets_generator(n_elem_layer_list):
    # n_elem_layer_list = [2,3,4] will produce all 2-aspect nets with
    # all nodelayers {0,1} x {0,1,2} x {0,1,2,3} and all edge configurations
    # (i.e. all nodelayers are in and powerset of possible edges).
    M_full = pymnet.MultilayerNetwork(aspects=len(n_elem_layer_list)-1,fullyInterconnected=False)
    nl_list = list(itertools.product(*[range(n) for n in n_elem_layer_list]))
    for ii,nl1 in enumerate(nl_list):
        for nl2 in nl_list[ii+1:]:
            M_full[nl1][nl2] = 1
    for m in pymnet.transforms.subnet_iter(M_full):
        yield m

def read_temp_file_res(filename):
    # returns tuple ((nl-mesu-runtime,nl-mesu-n-subnets),(a-mesu-runtime,a-mesu-n-subnets))
    # order defined by order of writing in result file
    res = []
    with open(filename,'r') as f:
        for line in f:
            line_arr = line.strip().split()
            curr_record = (float(line_arr[1]),int(line_arr[2]))
            res.append(curr_record)
    return tuple(res)

def call_temp_file(M,subnet_sizes):
    temp_dir = '/tmp/'
    inputfile = temp_dir + 'net_temp.edges'
    save_edgelist_cpp_format(M, inputfile)
    for subnet_size in subnet_sizes:
        n_aspects = len(subnet_size)-1
        outputfile = temp_dir + 'results_temp'
        call_str = './mesu_' + str(n_aspects) + '.out' + ' ' + "'" + inputfile + "'" + ' ' + "'" + outputfile + "'" + ' ' + "'" + str(subnet_size).replace(' ','').strip('()') + "'"
        os.system(call_str)
        res = read_temp_file_res(outputfile)
        os.remove(outputfile)
        yield res,subnet_size
    os.remove(inputfile)

def run_exhaustive_search(n_elem_layer_list,result_dir,save_every=10**5):
    res_list = []
    anomaly_list = []
    # subnet_sizes everything where each size is at least 2
    # subnet_sizes = list(itertools.product(*[range(2,n+1) for n in n_elem_layer_list]))
    # too many results, run first with only max subnet_size?
    subnet_sizes = [tuple(n_elem_layer_list)]
    current_block_start = 0
    for fname in os.listdir(result_dir):
        fname_list = fname.split('_')
        current_block_start = max(current_block_start,int(fname_list[1])+1)
    current_block_end = current_block_start + save_every - 1
    ii = current_block_start
    init_block_start = current_block_start
    # skip to current block start
    for m in itertools.islice(make_all_small_nets_generator(n_elem_layer_list),init_block_start,None):
        for res,subnet_size in call_temp_file(m,subnet_sizes):
            res_list.append((res[0][1],res[1][1]))
            if res[0][1] != res[1][1]:
                anomaly_list.append((m,subnet_size,res))
        if ii == current_block_end:
            with open(result_dir+str(current_block_start)+'_'+str(current_block_end),'wb') as f:
                pickle.dump(res_list,f)
            if anomaly_list:
                with open(result_dir+str(current_block_start)+'_'+str(current_block_end)+'_anomaly','wb') as g:
                    pickle.dump(anomaly_list,g)
            current_block_start = current_block_end + 1
            current_block_end = current_block_end + save_every
            res_list = []
            anomaly_list = []
        ii = ii+1
    with open(result_dir+str(current_block_start)+'_-1','wb') as h:
        pickle.dump(res_list,h)
    if anomaly_list:
        with open(result_dir+str(current_block_start)+'_-1_anomaly','wb') as j:
            pickle.dump(anomaly_list,j)

def summarize_results(result_dir,only_binary=False):
    problem_list = []
    for fname in sorted(os.listdir(result_dir)):
        if 'anomaly' in fname:
            print('Anomaly found (file ' + fname + ')')
        else:
            with open(result_dir + fname,'rb') as f:
                res_list = pickle.load(f)
            for r in res_list:
                if r[0] != r[1]:
                    problem_list.append(fname)
                if only_binary and not (r[0] == 0 or r[0] == 1):
                    print('Value other than 0 or 1 encountered at ' + fname)
                    problem_list.append(fname)
    if problem_list:
        print('Problems encountered, see return value')
    else:
        print('No problems encountered')
    return problem_list

if __name__ == '__main__':
    n_elem_layer_list = [int(x) for x in sys.argv[1].strip().split(',')]
    result_dir = sys.argv[2]
    if len(sys.argv) > 3:
        save_every = int(sys.argv[3])
    else:
        save_every = 10**5
    run_exhaustive_search(n_elem_layer_list,result_dir,save_every)
