#include <iostream>                  // for std::cout
#include <fstream>                   // for file reading
#include <filesystem>                // for checking if file exists
#include <string>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <array>                     // for copyable elementary layer arrays
#include <typeinfo>
#include <functional>                // hash template
#include <boost/container_hash/hash.hpp> // more hash functions from boost
#include <unordered_map>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <chrono>                    // for timing
#define N_ASPECTS 2

using namespace boost;

// multilayer data structures ------------------------------------------------------------------------------------------------

// node-layer class
class NL {
    std::array<int,N_ASPECTS+1> elem_layers;
    public:
     NL() {}; // default constructor
     NL(std::array<int,N_ASPECTS+1> el) : elem_layers(el) {}
     std::array<int,N_ASPECTS+1> get_el(void) const {return elem_layers;}
     // for testing shallow vs deep copy (i.e. is the array deepcopied)
     void set_el(std::array<int,N_ASPECTS+1> new_el) {elem_layers = new_el;}
     void print(void) const {std::cout << "{"; for (int ii=0; ii<elem_layers.size(); ii++) {std::cout << elem_layers[ii]; if (ii != elem_layers.size()-1) std::cout << ",";} std::cout << "}";}
     bool operator==(const NL& other) const {return elem_layers == other.elem_layers;}
     // are explicit copy and move constructors required for unordered_map?
};

// hash for NL (hash elem_layers)
namespace std {
    template <> struct hash<NL> {
        std::size_t operator()(const NL& nl) const {
            return hash_value(nl.get_el());
        }
    };
}

using Graph = adjacency_list<setS, vecS, undirectedS>;
using Vertex =  graph_traits<Graph>::vertex_descriptor;
using VertexIterator = graph_traits<Graph>::vertex_iterator;
using EdgeIterator = graph_traits<Graph>::edge_iterator;
using AdjacencyIterator = graph_traits<Graph>::adjacency_iterator;
using Degree = graph_traits<Graph>::degree_size_type;

class MLnet {
    Graph m;
    // Key, value
    std::unordered_map<NL,Vertex> nls;
    std::unordered_map<Vertex,NL> reverse_nls;
    public:
     int count_nodelayer(NL nl) const {return nls.count(nl);}
     void add_nodelayer(NL nl) {
        if (nls.count(nl) < 1) {
            Vertex vertex_descriptor = add_vertex(m);
            nls[nl] = vertex_descriptor;
            reverse_nls[vertex_descriptor] = nl;
        }
     }
     void add_nodelayer(std::array<int,N_ASPECTS+1> elem_layers) {
        NL nl (elem_layers);
        add_nodelayer(nl);
     }
     Vertex get_id_from_nl(NL nl) const {return nls.at(nl);}
     Vertex get_id_from_nl(std::array<int,N_ASPECTS+1> elem_layers) const {NL nl (elem_layers);return nls.at(nl);}
     NL get_nl_from_id(Vertex vertex) const {return reverse_nls.at(vertex);}
     // iterate through unordered_map nls to find all nodelayers (i.e. not implemented from BGL Graph)
     std::pair<std::vector<NL>,std::vector<Vertex>> get_all_nls() const {
        std::vector<NL> nodelayers;
        nodelayers.reserve(nls.size());
        std::vector<Vertex> vertices;
        vertices.reserve(nls.size());
        for(auto kv : nls) {
            nodelayers.push_back(kv.first);
            vertices.push_back(kv.second);
        }
        std::pair<std::vector<NL>,std::vector<Vertex>> combined;
        combined = std::make_pair(nodelayers,vertices);
        return combined;
     }
     void print_all_nls() const {
        std::pair<std::vector<NL>,std::vector<Vertex>> combined = get_all_nls();
        for (int i=0; i<combined.first.size(); i++) {
            combined.first[i].print();
            std::cout << ", id: " << combined.second[i];
            std::cout << "\n";
        }
     }
     // get all nls using the Boost underlying graph; for testing purposes
     void print_all_nls_from_underlying_graph() const {
        std::pair<VertexIterator,VertexIterator> all_nl_ids = vertices(m);
        for (Vertex nl_id : make_iterator_range(all_nl_ids)) {
            std::cout << nl_id << "\n";
        }
     }
     // has to be not named add_edge, because scope resolution happens before overload resolution?
     void add_mledge(NL nl1, NL nl2) {
        Vertex vert1 = get_id_from_nl(nl1);
        Vertex vert2 = get_id_from_nl(nl2);
        add_edge(vert1,vert2,m);
     }
     void add_mledge(std::array<int,N_ASPECTS+1> elem_layers1, std::array<int,N_ASPECTS+1> elem_layers2) {
        NL nl1 (elem_layers1);
        NL nl2 (elem_layers2);
        add_mledge(nl1,nl2);
     }
     // fill with all possible edges (not self-edges)
     void fill_mledges() {
        std::pair<std::vector<NL>,std::vector<Vertex>> combined = get_all_nls();
        for (int ii=0; ii<combined.first.size(); ii++) {
            for (int jj=ii+1; jj<combined.first.size(); jj++) {add_mledge(combined.first[ii],combined.first[jj]);}
        }
     }
     std::pair<EdgeIterator,EdgeIterator> get_all_mledges() const {return edges(m);}
     void print_all_mledges() const {
        std::pair<EdgeIterator,EdgeIterator> edge_iter_pair = get_all_mledges();
        for (EdgeIterator edge = edge_iter_pair.first; edge != edge_iter_pair.second; edge++) {
            Vertex src = source(*edge,m);
            Vertex trg = target(*edge,m);
            NL src_nl = get_nl_from_id(src);
            NL trg_nl = get_nl_from_id(trg);
            std::cout << "("; src_nl.print(); std::cout << ", "; trg_nl.print(); std::cout << ")";
            std::cout << " ids: " << src << ", " << trg << "\n";
        }
     }
     std::vector<NL> get_neighbors(NL nl) const {
        Vertex id = get_id_from_nl(nl);
        std::vector<NL> neighbors;
        neighbors.reserve(degree(id,m));
        std::pair<AdjacencyIterator,AdjacencyIterator> adj = adjacent_vertices(id,m);
        for (Vertex neighbor : make_iterator_range(adj)) {
            neighbors.push_back(get_nl_from_id(neighbor));
        }
        return neighbors;
     }
     std::vector<NL> get_neighbors(std::array<int,N_ASPECTS+1> elem_layers) const {NL nl (elem_layers); return get_neighbors(nl);}
     void print_neighbors(NL nl) const {
        std::vector<NL> neighbors = get_neighbors(nl);
        for (NL neighbor : neighbors) {neighbor.print(); if (!(neighbor == neighbors.back())) std::cout << ", ";}
     }
     void print_neighbors(std::array<int,N_ASPECTS+1> elem_layers) const {
        NL nl (elem_layers);
        print_neighbors(nl);
     }
     Degree get_degree(NL nl) const {Vertex v = get_id_from_nl(nl); return degree(v,m);}
     Degree get_degree(std::array<int,N_ASPECTS+1> elem_layers) const {Vertex v = get_id_from_nl(elem_layers); return degree(v,m);}
     void print_all_degrees() const {
        std::pair<std::vector<NL>,std::vector<Vertex>> combined = get_all_nls();
        for (int i=0; i<combined.first.size(); i++) {
            combined.first[i].print();
            Degree deg = degree(combined.second[i],m);
            std::cout << ", degree: " << deg;
            std::cout << "\n";
        }
     }
     // subnet function: takes array of vectors of elementary layers as parameter
     // NB! Vertex IDs will probably not match original net (but nodelayers will)
     MLnet subnet(std::array<std::vector<int>,N_ASPECTS+1> subnet_elem_layers) const {
        MLnet new_subnet;
        // cartesian product of elementary layers in subnet_elem_layers to get all possible nodelayers
        // we get total length and then do divisions and modulos to get the proper index to iterate over all product elements
        int total_length = 1;
        std::array<int,N_ASPECTS+1> divisors;
        std::array<int,N_ASPECTS+1> modulos;
        for (int ii=0; ii<N_ASPECTS+1; ii++) {
            total_length = total_length * subnet_elem_layers[ii].size();
            modulos[ii] = subnet_elem_layers[ii].size();
            divisors[ii] = 1;
            for (int jj=ii+1; jj<N_ASPECTS+1; jj++) {divisors[ii] = divisors[ii]*subnet_elem_layers[jj].size();}
        }
        for (int ii=0; ii<total_length; ii++) {
            std::array<int,N_ASPECTS+1> current_nodelayer;
            for (int jj=0; jj<N_ASPECTS+1; jj++) {current_nodelayer[jj] = subnet_elem_layers[jj][(ii/divisors[jj])%modulos[jj]];}
            // finally we get the possible nodelayers in the subnet
            NL curr_nl = NL(current_nodelayer);
            if (count_nodelayer(curr_nl) > 0) {new_subnet.add_nodelayer(curr_nl);}
            //curr_nl.print(); std::cout << "\n";
        }
        // create edges
        // iterating over all neighbors
        // TODO: iterating over other nls if degree is high
        std::pair<std::vector<NL>,std::vector<Vertex>> combined = new_subnet.get_all_nls();
        for (int ii=0; ii<combined.first.size(); ii++) {
            std::vector<NL> neighbors = get_neighbors(combined.first[ii]);
            for (auto neigh : neighbors) {
                if (new_subnet.count_nodelayer(neigh) > 0) {new_subnet.add_mledge(combined.first[ii],neigh);}
            }
        }
        return new_subnet;
     }
     bool is_connected() const {
        std::vector<int> component(num_vertices(m));
        int number_of_components = connected_components(m, &component[0]);
        if (number_of_components == 1) {return true;}
        else {return false;}
     }
};

// algorithms ----------------------------------------------------------------------------------------------------------------

// NL-MESU

std::array<std::unordered_set<int>,N_ASPECTS+1> spanned_space(const std::unordered_set<NL>& VM) {
    std::array<std::unordered_set<int>,N_ASPECTS+1> space;
    for (auto vm_iter = VM.begin(); vm_iter != VM.end(); vm_iter++) {
        std::array<int,N_ASPECTS+1> elem_layers = (*vm_iter).get_el();
        for (int ii=0; ii<elem_layers.size(); ii++) {
            space[ii].insert(elem_layers[ii]);
        }
    }
    return space;
}

std::array<int,N_ASPECTS+1> spanned_volume(const std::unordered_set<NL>& VM) {
    std::array<std::unordered_set<int>,N_ASPECTS+1> space = spanned_space(VM);
    std::array<int,N_ASPECTS+1> volume;
    for (int jj=0; jj<space.size(); jj++) {
        volume[jj] = space[jj].size();
    }
    return volume;
}

std::unordered_set<NL> VM_neighbors(const MLnet& mlnet, const std::unordered_set<NL>& VM) {
    std::unordered_set<NL> all_neighs;
    for (auto vm_iter = VM.begin(); vm_iter != VM.end(); vm_iter++) {
        std::vector<NL> curr_neighs = mlnet.get_neighbors(*vm_iter);
        for (NL neighbor : curr_neighs) {
            if (VM.count(neighbor) < 1) {all_neighs.insert(neighbor);}
        }
    }
    return all_neighs;
}

bool valid_nl_mesu(const MLnet& mlnet, const std::unordered_set<NL>& VM_subnet, const std::unordered_set<NL>& extension, const Vertex& gamma_index) {
    std::array<std::unordered_set<int>,N_ASPECTS+1> space = spanned_space(VM_subnet);
    // change from set into vector so subnet index magic works
    std::array<std::vector<int>,N_ASPECTS+1> subnet_elem_layers;
    for (int ii=0; ii<N_ASPECTS+1; ii++) {
        subnet_elem_layers[ii].reserve(space[ii].size());
        for (auto it=space[ii].begin(); it!=space[ii].end(); ) {
            subnet_elem_layers[ii].push_back(std::move(space[ii].extract(it++).value())); // this modifies space -> iterator invalid?
        }
    }
    MLnet sub = mlnet.subnet(subnet_elem_layers);
    if (sub.is_connected()) {
        bool all_indices_valid = true;
        std::pair<std::vector<NL>,std::vector<Vertex>> all_nls = sub.get_all_nls();
        // get id from original net, not subnet
        for (NL sub_nl : all_nls.first) {if (mlnet.get_id_from_nl(sub_nl) < gamma_index) {all_indices_valid = false;break;}}
        if (all_indices_valid) {
            std::unordered_set<NL> VM_neighs = VM_neighbors(sub,VM_subnet);
            bool all_neighs_in_extension = true;
            for (NL neigh : VM_neighs) {if (extension.count(neigh) < 1) {all_neighs_in_extension = false;break;}}
            if (all_neighs_in_extension) {
                return true;
            }
        }
    }
    return false;
}

void extend_nl_mesu(const MLnet& mlnet, const std::array<int,N_ASPECTS+1> size, std::unordered_set<NL>& VM_subnet, std::unordered_set<NL>& extension, Vertex& gamma_index, int& total_number) {
    std::array<int,N_ASPECTS+1> volume = spanned_volume(VM_subnet);
    if (volume == size) {
        if (valid_nl_mesu(mlnet,VM_subnet,extension,gamma_index)) {total_number++;}
        return;
    }
    std::unordered_set<NL> N = VM_neighbors(mlnet,VM_subnet);
    while (not extension.empty()) {
        // pop element; no in-built popping function available
        NL gamma_prime = *extension.begin();
        extension.erase(extension.begin());
        // new VM_subnet
        std::unordered_set<NL> VM_subnet_prime = VM_subnet;
        VM_subnet_prime.insert(gamma_prime);
        std::array<int,N_ASPECTS+1> subnet_prime_volume = spanned_volume(VM_subnet_prime);
        // check if any aspect goes over the size limit; no easy solution so just use a dumb flag approach
        bool oversize = false;
        for (int ii=0; ii<size.size(); ii++) {if (subnet_prime_volume[ii]>size[ii]) {oversize=true;break;}}
        if (oversize) {continue;}
        std::unordered_set<NL> extension_prime = extension;
        std::vector<NL> neighbors = mlnet.get_neighbors(gamma_prime);
        for (NL neigh : neighbors) {
            if (VM_subnet.count(neigh) < 1 and N.count(neigh) < 1 and mlnet.get_id_from_nl(neigh) > gamma_index) {extension_prime.insert(neigh);}
        }
        extend_nl_mesu(mlnet,size,VM_subnet_prime,extension_prime,gamma_index,total_number);
    }
    return;
}

int nl_mesu(const MLnet& mlnet, const std::array<int,N_ASPECTS+1> size) {
    int total_number = 0;
    std::pair<std::vector<NL>,std::vector<Vertex>> combined = mlnet.get_all_nls();
    for (int ii=0; ii<combined.first.size(); ii++) {
        Vertex gamma_index = mlnet.get_id_from_nl(combined.first[ii]);
        std::unordered_set<NL> extension;
        std::unordered_set<NL> VM_subnet;
        VM_subnet.insert(combined.first[ii]);
        std::vector<NL> neighbors = mlnet.get_neighbors(combined.first[ii]);
        for (NL neigh : neighbors) {
            if (mlnet.get_id_from_nl(neigh) > gamma_index) {
                extension.insert(neigh);
            }
        }
        //std::cout << "nl: ";
        //combined.first[ii].print();
        //std::cout << " ext: ";
        //for (auto e = extension.begin(); e != extension.end(); e++) {(*e).print();}
        //std::cout << "\n";
        extend_nl_mesu(mlnet,size,VM_subnet,extension,gamma_index,total_number);
    }
    return total_number;
}

// A-MESU

// cartesian product copied from subnet, casting into vector copied from valid_nl_mesu; maybe a better way to not duplicate code?
std::unordered_set<NL> subnet_nodelayers(const MLnet& mlnet, const std::array<std::unordered_set<int>,N_ASPECTS+1>& S) {
    std::unordered_set<NL> sub_nls;
    int total_length = 1;
    std::array<int,N_ASPECTS+1> divisors;
    std::array<int,N_ASPECTS+1> modulos;
    std::array<std::vector<int>,N_ASPECTS+1> subnet_elem_layers;
    // ii and jj of subnet_elem_layers required, need to loop twice
    for (int ii=0; ii<N_ASPECTS+1; ii++) {
        subnet_elem_layers[ii].reserve(S[ii].size());
        for (auto it=S[ii].begin(); it!=S[ii].end(); it++) {
            subnet_elem_layers[ii].push_back(*it);
        }
    }
    for (int jj=0; jj<N_ASPECTS+1; jj++) {
        total_length = total_length * subnet_elem_layers[jj].size();
        modulos[jj] = subnet_elem_layers[jj].size();
        divisors[jj] = 1;
        for (int kk=jj+1; kk<N_ASPECTS+1; kk++) {divisors[jj] = divisors[jj]*subnet_elem_layers[kk].size();}
    }
    for (int ll=0; ll<total_length; ll++) {
        std::array<int,N_ASPECTS+1> current_nodelayer;
        for (int mm=0; mm<N_ASPECTS+1; mm++) {current_nodelayer[mm] = subnet_elem_layers[mm][(ll/divisors[mm])%modulos[mm]];}
        // finally we get the possible nodelayers in the subnet
        NL curr_nl = NL(current_nodelayer);
        if (mlnet.count_nodelayer(curr_nl) > 0) {sub_nls.insert(curr_nl);}
    }
    return sub_nls;
}

std::unordered_set<NL> subnet_diff(const MLnet& mlnet, const std::array<std::unordered_set<int>,N_ASPECTS+1>& S_prime, const std::array<std::unordered_set<int>,N_ASPECTS+1>& S) {
    std::unordered_set<NL> nodelayers_S_prime = subnet_nodelayers(mlnet,S_prime);
    std::unordered_set<NL> nodelayers_S = subnet_nodelayers(mlnet,S);
    // remove in-place from nodelayers_S_prime
    for (NL nlS : nodelayers_S) {nodelayers_S_prime.erase(nlS);}
    return nodelayers_S_prime;
}

std::unordered_set<NL> subnet_valid_neighbor_nls(const MLnet& mlnet, const std::array<std::unordered_set<int>,N_ASPECTS+1>& S, const Vertex& gamma_index) {
    // cast sets into vectors
    std::unordered_set<NL> sub_nls = subnet_nodelayers(mlnet, S);
    std::unordered_set<NL> valid_neighs;
    for (NL nl : sub_nls) {
        for (NL neigh : mlnet.get_neighbors(nl)) {
            if (sub_nls.count(neigh) < 1 and valid_neighs.count(neigh) < 1) {
                bool neigh_indices_valid = true;
                std::array<std::unordered_set<int>,N_ASPECTS+1> S_prime = S;
                std::array<int,N_ASPECTS+1> neigh_el = neigh.get_el();
                for (int jj=0; jj<N_ASPECTS+1; jj++) {S_prime[jj].insert(neigh_el[jj]);}
                for (NL possible_addition : subnet_diff(mlnet,S_prime,S)) {
                    if (mlnet.get_id_from_nl(possible_addition) < gamma_index) {
                        neigh_indices_valid = false;
                        break;
                    }
                }
                if (neigh_indices_valid) {valid_neighs.insert(neigh);}
            }
        }
    }
    return valid_neighs;
}

std::array<std::unordered_set<int>,N_ASPECTS+1> subnet_valid_neighbor_elem_layers(const MLnet& mlnet, const std::array<std::unordered_set<int>,N_ASPECTS+1>& S, const Vertex& gamma_index) {
    std::array<std::unordered_set<int>,N_ASPECTS+1> N;
    for (NL valid_neigh : subnet_valid_neighbor_nls(mlnet,S,gamma_index)) {
        std::array<int,N_ASPECTS+1> neigh_el = valid_neigh.get_el();
        for (int ii=0; ii<N_ASPECTS+1; ii++) {N[ii].insert(neigh_el[ii]);}
    }
    return N;
}

bool valid_a_mesu(const MLnet& mlnet, const std::array<std::unordered_set<int>,N_ASPECTS+1>& S) {
    // cast sets into vectors
    std::array<std::vector<int>,N_ASPECTS+1> subnet_elem_layers;
    for (int ii=0; ii<N_ASPECTS+1; ii++) {
        subnet_elem_layers[ii].reserve(S[ii].size());
        for (int el : S[ii]) {subnet_elem_layers[ii].push_back(el);}
    }
    MLnet sub = mlnet.subnet(subnet_elem_layers);
    if (sub.is_connected()) {
        // iterate over nl's and make valid_S from their elem layers, compare to S
        std::pair<std::vector<NL>,std::vector<Vertex>> all_nls = sub.get_all_nls();
        std::array<std::unordered_set<int>,N_ASPECTS+1> S_valid;
        for (NL nl : all_nls.first) {std::array<int,N_ASPECTS+1> el=nl.get_el(); for (int ii=0; ii<N_ASPECTS+1; ii++) {S_valid[ii].insert(el[ii]);}}
        if (S == S_valid) {return true;}
    }
    return false;
}

std::vector<int> candidate_extension_indices(const std::array<std::unordered_set<int>,N_ASPECTS+1>& extension, std::array<std::unordered_set<int>,N_ASPECTS+1>& S, const std::array<int,N_ASPECTS+1> size) {
    std::vector<int> candidates;
    for (int ii=0; ii<N_ASPECTS+1; ii++) {if ((not extension[ii].empty()) and (S[ii].size() < size[ii])) {candidates.push_back(ii);}}
    return candidates;
}

void extend_a_mesu(const MLnet& mlnet, const std::array<int,N_ASPECTS+1> size, std::array<std::unordered_set<int>,N_ASPECTS+1>& S, std::array<std::unordered_set<int>,N_ASPECTS+1>& extension, Vertex& gamma_index, int& total_number) {
    bool size_ok = true;
    for (int ii=0; ii<N_ASPECTS+1; ii++) {if (S[ii].size()!=size[ii]) {size_ok=false;break;}}
    if (size_ok) {
        if (valid_a_mesu(mlnet, S)) {total_number++;}
        return;
    }
    std::array<std::unordered_set<int>,N_ASPECTS+1> N = subnet_valid_neighbor_elem_layers(mlnet,S,gamma_index);
    std::vector<int> possible_indices = candidate_extension_indices(extension,S,size);
    while (not possible_indices.empty()) {
        int chosen_index = possible_indices.front();
        // get l
        int l = *(extension[chosen_index].begin());
        extension[chosen_index].erase(extension[chosen_index].begin());
        std::array<std::unordered_set<int>,N_ASPECTS+1> S_prime;
        std::array<std::unordered_set<int>,N_ASPECTS+1> extension_prime;
        // is this deepcopy?
        for (int ii=0; ii<N_ASPECTS+1; ii++) {S_prime[ii] = S[ii]; extension_prime[ii] = extension[ii];}
        S_prime[chosen_index].insert(l);
        std::unordered_set<NL> S_prime_nls = subnet_nodelayers(mlnet,S_prime);
        // make new possible_indices
        possible_indices = candidate_extension_indices(extension,S,size);
        for (NL tau : subnet_diff(mlnet,S_prime,S)) {
            for (NL delta : mlnet.get_neighbors(tau)) {
                if (S_prime_nls.count(delta) < 1) {
                    std::array<std::unordered_set<int>,N_ASPECTS+1> S_prime_with_delta = S_prime;
                    std::array<int,N_ASPECTS+1> delta_els = delta.get_el();
                    for (int jj=0; jj<N_ASPECTS+1; jj++) {S_prime_with_delta[jj].insert(delta_els[jj]);}
                    bool lambda_indices_valid = true;
                    for (NL lambda : subnet_diff(mlnet,S_prime_with_delta,S_prime)) {
                        if (mlnet.get_id_from_nl(lambda) < gamma_index) {lambda_indices_valid = false;break;}
                    }
                    if (lambda_indices_valid) {for (int kk=0; kk<N_ASPECTS+1; kk++) {if (N[kk].count(delta_els[kk]) < 1 and S_prime[kk].count(delta_els[kk]) < 1) {extension_prime[kk].insert(delta_els[kk]);}}}
                }
            }
        }
        extend_a_mesu(mlnet,size,S_prime,extension_prime,gamma_index,total_number);
    }
}

int a_mesu(const MLnet& mlnet, const std::array<int,N_ASPECTS+1> size) {
    int total_number = 0;
    std::pair<std::vector<NL>,std::vector<Vertex>> combined = mlnet.get_all_nls();
    for (int ii=0; ii<combined.first.size(); ii++) {
        Vertex gamma_index = mlnet.get_id_from_nl(combined.first[ii]);
        std::array<std::unordered_set<int>,N_ASPECTS+1> S;
        for (int jj=0; jj<N_ASPECTS+1; jj++) {S[jj].insert(combined.first[ii].get_el()[jj]);}
        std::array<std::unordered_set<int>,N_ASPECTS+1> extension;
        std::vector<NL> neighbors = mlnet.get_neighbors(combined.first[ii]);
        for (NL neigh : neighbors) {
            std::array<std::unordered_set<int>,N_ASPECTS+1> S_prime = S;
            for (int jj=0; jj<N_ASPECTS+1; jj++) {S_prime[jj].insert(neigh.get_el()[jj]);}
            std::unordered_set<NL> sub_diff = subnet_diff(mlnet,S_prime,S);
            bool lambda_indices_valid = true;
            for (NL lambda : sub_diff) {if (mlnet.get_id_from_nl(lambda) < gamma_index) {lambda_indices_valid=false;break;}}
            if (lambda_indices_valid) {for (int jj=0; jj<N_ASPECTS+1; jj++) {int el=neigh.get_el()[jj]; if (S[jj].count(el)<1) {extension[jj].insert(el);}}}
        }
        extend_a_mesu(mlnet,size,S,extension,gamma_index,total_number);
    }
    return total_number;
}

// file io -------------------------------------------------------------------------------------------------------------------

MLnet load_edge_file(const std::string& filename) {
    // file format: a0 a1 a2 ... b0 b1 b2 ... (anything after ignored)
    // elementary layers need to be integers
    // only specifying one nodelayer per line adds that nodelayer
    MLnet mlnet;
    std::ifstream file(filename);
    std::string current_line;
    std::string current_value_str;
    std::array<std::array<int,N_ASPECTS+1>,2> current_elem_layers;
    int aspect_counter = 0;
    int edge_endpoint_counter = 0;
    while (std::getline(file,current_line)) {
        auto sstream = std::istringstream(current_line);
        while (sstream >> current_value_str) {
            current_elem_layers[edge_endpoint_counter][aspect_counter] = std::stoi(current_value_str);
            aspect_counter++;
            if (aspect_counter > N_ASPECTS) {
                mlnet.add_nodelayer(current_elem_layers[edge_endpoint_counter]);
                aspect_counter = 0;
                edge_endpoint_counter++;
            }
            if (edge_endpoint_counter > 1) {
                mlnet.add_mledge(current_elem_layers[0],current_elem_layers[1]);
                break;
            }
        }
        aspect_counter = 0;
        edge_endpoint_counter = 0;
    }
    file.close();
    return mlnet;
}

void run_edge_file(const std::string& filename, const std::string& savename, const std::array<int,N_ASPECTS+1>& size) {
    if (not std::filesystem::exists(savename)) {
        MLnet mlnet = load_edge_file(filename);
        std::ofstream output(savename);
        auto start = std::chrono::steady_clock::now();
        int number_of_subnets = nl_mesu(mlnet,size);
        auto end = std::chrono::steady_clock::now();
        auto tdiff = end-start;
        output << "nl-mesu "<<std::chrono::duration<double> (tdiff).count() << " " << number_of_subnets << "\n";
        start = std::chrono::steady_clock::now();
        number_of_subnets = a_mesu(mlnet,size);
        end = std::chrono::steady_clock::now();
        tdiff = end-start;
        output << "a-mesu "<<std::chrono::duration<double> (tdiff).count() << " " << number_of_subnets << "\n";
        output.close();
    }
}

std::array<int,N_ASPECTS+1> parse_size(const std::string& size_str, const char& delimiter) {
    auto sstream = std::istringstream(size_str);
    std::array<int,N_ASPECTS+1> elem_layers;
    std::string current_elem_layer;
    int aspect_counter = 0;
    while (std::getline(sstream,current_elem_layer,delimiter)) {
        elem_layers[aspect_counter] = std::stoi(current_elem_layer);
        aspect_counter++;
    }
    assert(N_ASPECTS+1==aspect_counter);
    return elem_layers;
}

MLnet load_ppi_data(const std::string filename) {
    MLnet mlnet;
    std::ifstream file(filename);
    int layerID, node1ID, node2ID, edgeweight;
    std::unordered_set<int> nodes;
    std::set<int> layers;
    std::unordered_map<int,std::unordered_set<int>> nodes_by_layer;
    while (file >> layerID >> node1ID >> node2ID >> edgeweight) {
        mlnet.add_nodelayer({node1ID,layerID});
        mlnet.add_nodelayer({node2ID,layerID});
        mlnet.add_mledge({node1ID,layerID},{node2ID,layerID});
        nodes_by_layer[layerID].insert(node1ID);
        nodes_by_layer[layerID].insert(node2ID);
        nodes.insert(node1ID);
        nodes.insert(node2ID);
        layers.insert(layerID);
    }
    file.close();
    // add interlayer edges
    // (NB! only nodelayers with connections on a layer are included on that layer)
    for (int node : nodes) {
        for (auto it1 = layers.begin(); it1 != layers.end(); it1++) {
            for (auto it2 = layers.begin(); it2 != layers.end(); it2++) {
                if (*it1 < *it2 and nodes_by_layer[*it1].count(node) > 0 and nodes_by_layer[*it2].count(node) > 0) {mlnet.add_mledge({node,*it1},{node,*it2});}
            }
        }
    }
    return mlnet;
}

void run_ppi_data(const std::string& filename, const std::string& savename, const std::array<int,N_ASPECTS+1>& size) {
    if (not std::filesystem::exists(savename)) {
        MLnet mlnet = load_ppi_data(filename);
        std::ofstream output(savename);
        auto start = std::chrono::steady_clock::now();
        int number_of_subnets = nl_mesu(mlnet,size);
        auto end = std::chrono::steady_clock::now();
        auto tdiff = end-start;
        output << "nl-mesu "<<std::chrono::duration<double> (tdiff).count() << " " << number_of_subnets << "\n";
        start = std::chrono::steady_clock::now();
        number_of_subnets = a_mesu(mlnet,size);
        end = std::chrono::steady_clock::now();
        tdiff = end-start;
        output << "a-mesu "<<std::chrono::duration<double> (tdiff).count() << " " << number_of_subnets << "\n";
        output.close();
    }
}

void run_all_ppi(int specific_index = -1) {
    std::vector<std::string> filenames = {"multiplex_pp_data/Arabidopsis_Multiplex_Genetic/Dataset/arabidopsis_genetic_multiplex.edges",
                                          "multiplex_pp_data/Bos_Multiplex_Genetic/Dataset/bos_genetic_multiplex.edges",
                                          "multiplex_pp_data/Candida_Multiplex_Genetic/Dataset/candida_genetic_multiplex.edges",
                                          "multiplex_pp_data/Celegans_Multiplex_Genetic/Dataset/celegans_genetic_multiplex.edges",
                                          "multiplex_pp_data/Drosophila_Multiplex_Genetic/Dataset/drosophila_genetic_multiplex.edges",
                                          "multiplex_pp_data/Gallus_Multiplex_Genetic/Dataset/gallus_genetic_multiplex.edges",
                                          "multiplex_pp_data/Mus_Multiplex_Genetic/Dataset/mus_genetic_multiplex.edges",
                                          "multiplex_pp_data/Plasmodium_Multiplex_Genetic/Dataset/plasmodium_genetic_multiplex.edges",
                                          "multiplex_pp_data/Rattus_Multiplex_Genetic/Dataset/rattus_genetic_multiplex.edges",
                                          "multiplex_pp_data/SacchCere_Multiplex_Genetic/Dataset/sacchcere_genetic_multiplex.edges",
                                          "multiplex_pp_data/SacchPomb_Multiplex_Genetic/Dataset/sacchpomb_genetic_multiplex.edges"
                                         };
    std::vector<std::string> ids = {"arabidopsis",
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
                                   };
    std::vector<std::array<int,N_ASPECTS+1>> sizes = {{2,2},{3,2},{2,3},{3,3},{4,2},{4,3}};
    std::string savename;
    for (std::array<int,N_ASPECTS+1> size : sizes) {
        for (int ii=0; ii<filenames.size(); ii++) {
            savename = "cpp_results/" + ids[ii] + "_(" + std::to_string(size[0]) + "," + std::to_string(size[1]) + ").txt";
            if (specific_index == -1 or specific_index == ii) {
                run_ppi_data(filenames[ii],savename,size);
                std::cout << ids[ii] << " (" << size[0] << "," << size[1] << ") done" << std::endl;
            }
        }
    }
}

void print_arguments(const std::vector<std::string>& args) {
    for (std::string param_str : args) {
        std::cout << param_str << "\n";
    }
}

void test_write_to_argument_file(const std::string& fname) {
    std::ofstream output(fname);
    output << "Writing text to file\n";
    output << "and on the second line\n";
    output.close();
}

// main ----------------------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    /*
    MLnet mlnet;
    mlnet.add_nodelayer({5,6,7});
    mlnet.add_nodelayer({2,3,4});
    mlnet.add_nodelayer({1,2,3});
    std::cout << "Nodelayers in net: \n";
    mlnet.print_all_nls();
    std::cout << "Testing Boost vertices function to see that the nodelayers are actually there: \n";
    mlnet.print_all_nls_from_underlying_graph();
    mlnet.add_mledge({2,3,4},{5,6,7});
    mlnet.add_mledge({2,3,4},{1,2,3});
    std::cout << "Edges:" << "\n";
    mlnet.print_all_mledges();
    std::cout << "Testing access to neighbors of {2,3,4} from curly bracket notation: ";
    auto nn = mlnet.get_neighbors({2,3,4});
    for (auto nx : nn) nx.print(); std::cout << "\n";
    std::cout << "Neighbors of {2,3,4}: ";
    mlnet.print_neighbors({2,3,4});
    std::cout << "\n";
    std::cout << "Neighbors of {5,6,7}: ";
    mlnet.print_neighbors({5,6,7});
    std::cout << "\n";
    Degree d = mlnet.get_degree({2,3,4});
    std::cout << "Degree of {2,3,4}: " << d << "\n";
    std::cout << "All degrees:\n";
    mlnet.print_all_degrees();
    std::array<std::vector<int>,N_ASPECTS+1> subnet_els;
    subnet_els[0] = std::vector<int> {1,2,3};
    subnet_els[1] = std::vector<int> {2,3};
    subnet_els[2] = std::vector<int> {3,4,5,6};
    MLnet sub = mlnet.subnet(subnet_els);
    std::cout << "Subnet nls:\n";
    sub.print_all_nls();
    std::cout << "Subnet edges:\n";
    sub.print_all_mledges();
    std::cout << "Subnet neighbors of {2,3,4}: ";
    sub.print_neighbors({2,3,4});
    std::cout << "\n";
    std::cout << "Subnet neighbors of {1,2,3}: ";
    sub.print_neighbors({1,2,3});
    std::cout << "\n";
    std::cout << "Connectedness:\n";
    if (mlnet.is_connected()) {std::cout << "mlnet is connected\n";} else {std::cout << "mlnet is not connected\n";}
    if (sub.is_connected()) {std::cout << "subnet is connected\n";} else {std::cout << "subnet is not connected\n";}
    std::cout << "Add isolated nl {3,2,1} to mlnet:\n";
    mlnet.add_nodelayer({3,2,1});
    if (mlnet.is_connected()) {std::cout << "mlnet is connected\n";} else {std::cout << "mlnet is not connected\n";}
    if (sub.is_connected()) {std::cout << "subnet is connected\n";} else {std::cout << "subnet is not connected\n";}
    std::cout << "Testing nl-mesu and a-mesu:\n";
    int number_of_subnets = nl_mesu(mlnet,{3,3,3});
    std::cout << "nl-mesu number of subnets: " << number_of_subnets << "\n";
    int number_of_subnets_2 = a_mesu(mlnet,{3,3,3});
    std::cout << "a-mesu number of subnets: " << number_of_subnets_2 << "\n";
    std::cout << "Add edge {5,6,7},{3,2,1}\n";
    mlnet.add_mledge({5,6,7},{3,2,1});
    number_of_subnets = nl_mesu(mlnet,{3,3,3});
    std::cout << "nl-mesu number of subnets: " << number_of_subnets << "\n";
    number_of_subnets_2 = a_mesu(mlnet,{3,3,3});
    std::cout << "a-mesu number of subnets: " << number_of_subnets_2 << "\n";
    std::cout << "Make a 2x2x2 complete network:\n";
    MLnet full;
    full.add_nodelayer({1,1,1});full.add_nodelayer({2,1,1});full.add_nodelayer({1,2,1});full.add_nodelayer({2,2,1});
    full.add_nodelayer({1,1,2});full.add_nodelayer({2,1,2});full.add_nodelayer({1,2,2});full.add_nodelayer({2,2,2});
    full.fill_mledges();
    std::cout << "Check nl-mesu and a-mesu:\n";
    number_of_subnets = nl_mesu(full,{2,2,2});
    number_of_subnets_2 = a_mesu(full,{2,2,2});
    std::cout << "Number of subnets {2,2,2}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 1)\n";
    number_of_subnets = nl_mesu(full,{2,1,2});
    number_of_subnets_2 = a_mesu(full,{2,1,2});
    std::cout << "Number of subnets {2,1,2}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 2)\n";
    number_of_subnets = nl_mesu(full,{1,1,2});
    number_of_subnets_2 = a_mesu(full,{1,1,2});
    std::cout << "Number of subnets {1,1,2}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 4)\n";
    number_of_subnets = nl_mesu(full,{1,1,1});
    number_of_subnets_2 = a_mesu(full,{1,1,1});
    std::cout << "Number of subnets {1,1,1}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 8)\n";
    number_of_subnets = nl_mesu(full,{1,1,3});
    number_of_subnets_2 = a_mesu(full,{1,1,3});
    std::cout << "Number of subnets {1,1,3}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 0)\n";
    std::cout << "Add nl {1,1,3} connected to {1,1,1}:\n";
    full.add_nodelayer({1,1,3});
    full.add_mledge({1,1,1},{1,1,3});
    number_of_subnets = nl_mesu(full,{1,2,2});
    number_of_subnets_2 = a_mesu(full,{1,2,2});
    std::cout << "Number of subnets {1,2,2}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 3)\n";
    number_of_subnets = nl_mesu(full,{2,1,3});
    number_of_subnets_2 = a_mesu(full,{2,1,3});
    std::cout << "Number of subnets {2,1,3}: nl: " << number_of_subnets << " a: " << number_of_subnets_2 << " (this number should be 1)\n";
    */
    /*
    auto start = std::chrono::steady_clock::now();
    MLnet celegans = load_ppi_data("multiplex_pp_data/Celegans_Multiplex_Genetic/Dataset/celegans_genetic_multiplex.edges");
    auto end = std::chrono::steady_clock::now();
    auto tdiff = end-start;
    std::cout << std::chrono::duration<double> (tdiff).count() << " s for loading net\n";
    // algos
    int number_of_subnets, number_of_subnets2;
    start = std::chrono::steady_clock::now();
    number_of_subnets = nl_mesu(celegans,{3,2});
    end = std::chrono::steady_clock::now();
    tdiff = end-start;
    std::cout << std::chrono::duration<double> (tdiff).count() << " s for 3,2 nl-mesu (" << number_of_subnets << " subnets)\n";
    start = std::chrono::steady_clock::now();
    number_of_subnets2 = a_mesu(celegans,{3,2});
    end = std::chrono::steady_clock::now();
    tdiff = end-start;
    std::cout << std::chrono::duration<double> (tdiff).count() << " s for 3,2 a-mesu (" << number_of_subnets2 << " subnets)\n";
    */
    //run_ppi_data("multiplex_pp_data/Celegans_Multiplex_Genetic/Dataset/celegans_genetic_multiplex.edges","cpp_results/celegans_(3,2).txt",{3,2});
    //run_all_ppi();
    std::vector<std::string> args(argv, argv+argc);
    print_arguments(args);
    if (argc > 1) {
        test_write_to_argument_file(args[1]);
    }
    if (argc > 2) {
        std::array<int,N_ASPECTS+1> subnet_size = parse_size(args[2],',');
        for (int ii : subnet_size) {std::cout << ii << "\n";}
    }
    //MLnet mlnet = load_edge_file("aaanetworkfiletestwrite.edges");
    //mlnet.print_all_nls();
    //mlnet.print_all_mledges();
}


