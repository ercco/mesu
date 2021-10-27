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
     void print(void) const {std::cout << "{"; for (int i : elem_layers) {std::cout << i; if (i != elem_layers.back()) std::cout << ",";} std::cout << "}";}
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

using Graph = adjacency_list<vecS, vecS, undirectedS>;
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
};

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
}


