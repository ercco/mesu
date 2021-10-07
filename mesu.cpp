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
     void print(void) const {for (int i : elem_layers) std::cout << i;}
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
};

int main() {
    MLnet mlnet;
    std::array<int,N_ASPECTS+1> a1 = {1,2,3};
    NL nl1 ({1,2,3});
    NL nl2;
    nl1.print();
    std::cout << "\n";
    nl2 = nl1;
    nl1.set_el({5,6,7});
    nl1.print();
    std::cout << "\n";
    nl2.print();
    std::cout << "\n";
    mlnet.add_nodelayer(nl1);
    mlnet.add_nodelayer({2,3,4});
    //Vertex x = mlnet.get_id_from_nl({1,2,3});
    //std::cout << "{1,2,3}: " << x << "\n";
    // we edited nl1 in the above code
    Vertex x = mlnet.get_id_from_nl({5,6,7});
    std::cout << "{5,6,7}: " << x << "\n";
    Vertex y = mlnet.get_id_from_nl({2,3,4});
    std::cout << "{2,3,4}: " << y << "\n";
    //Vertex z = mlnet.get_id_from_nl({3,4,5});
    //std::cout << "{3,4,5}: " << z << "\n";
}


