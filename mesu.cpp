#include <iostream>                  // for std::cout
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
#define N_ASPECTS 2

using namespace boost;

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
        Vertex vertex_descriptor = add_vertex(m);
        nls[nl] = vertex_descriptor;
        reverse_nls[vertex_descriptor] = nl;
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
        // TODO: extend_a_mesu
    }
    return total_number;
}

// main ----------------------------------------------------------------------------------------------------------------------

int main() {
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
    std::cout << "Testing nl-mesu:\n";
    int number_of_subnets = nl_mesu(mlnet,{3,3,3});
    std::cout << "Number of subnets: " << number_of_subnets << "\n";
    std::cout << "Add edge {5,6,7},{3,2,1}\n";
    mlnet.add_mledge({5,6,7},{3,2,1});
    number_of_subnets = nl_mesu(mlnet,{3,3,3});
    std::cout << "Number of subnets: " << number_of_subnets << "\n";
    std::cout << "Make a 2x2x2 complete network:\n";
    MLnet full;
    full.add_nodelayer({1,1,1});full.add_nodelayer({2,1,1});full.add_nodelayer({1,2,1});full.add_nodelayer({2,2,1});
    full.add_nodelayer({1,1,2});full.add_nodelayer({2,1,2});full.add_nodelayer({1,2,2});full.add_nodelayer({2,2,2});
    full.fill_mledges();
    std::cout << "Check nl-mesu:\n";
    number_of_subnets = nl_mesu(full,{2,2,2});
    std::cout << "Number of subnets {2,2,2}: " << number_of_subnets << " (this number should be 1)\n";
    number_of_subnets = nl_mesu(full,{2,1,2});
    std::cout << "Number of subnets {2,1,2}: " << number_of_subnets << " (this number should be 2)\n";
    number_of_subnets = nl_mesu(full,{1,1,2});
    std::cout << "Number of subnets {1,1,2}: " << number_of_subnets << " (this number should be 4)\n";
    number_of_subnets = nl_mesu(full,{1,1,1});
    std::cout << "Number of subnets {1,1,1}: " << number_of_subnets << " (this number should be 8)\n";
    number_of_subnets = nl_mesu(full,{1,1,3});
    std::cout << "Number of subnets {1,1,3}: " << number_of_subnets << " (this number should be 0)\n";
    std::cout << "Add nl {1,1,3} connected to {1,1,1}:\n";
    full.add_nodelayer({1,1,3});
    full.add_mledge({1,1,1},{1,1,3});
    number_of_subnets = nl_mesu(full,{1,2,2});
    std::cout << "Number of subnets {1,2,2}: " << number_of_subnets << " (this number should be 3)\n";
    number_of_subnets = nl_mesu(full,{2,1,3});
    std::cout << "Number of subnets {2,1,3}: " << number_of_subnets << " (this number should be 1)\n";
}


