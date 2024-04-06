// GraphConstructor.h
#ifndef GraphConstructor_HPP
#define GraphConstructor_HPP

#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include "ReadWrite.hpp"

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/all.hpp>
#include <vector>
#include <set>
#include <seqan3/alphabet/all.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits> // for std::decay_t
#include <boost/functional/hash.hpp>
#include <unordered_set>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graphml.hpp>
#include <seqan3/core/debug_stream.hpp> // for debug_stream
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/depth_first_search.hpp>

using namespace std;
using namespace seqan3::literals;
using namespace boost;

// Define the graph type with properties for vertices and edges
struct VertexProperties {
    std::vector<seqan3::dna5> read;
    uint32_t count;
};

struct EdgeProperties {
    // std::vector<seqan3::dna5> read1;
    // std::vector<seqan3::dna5> read2;    
    int weight;
};

namespace std {
template <>
struct hash<std::vector<seqan3::dna5>>
{
    std::size_t operator()(const std::vector<seqan3::dna5> & v) const
    {
        // Your custom hash implementation here, for example:
        std::size_t hash_value = 0;
        for (seqan3::dna5 character : v)
        {
            hash_value ^= seqan3::to_rank(character) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
        }
        return hash_value;
    }
};
}

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperties, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

template <class Name>
class custom_edge_label_writer {
public:
    custom_edge_label_writer(Name _name) : name(_name) {}
    template <class Edge>
    void operator()(std::ostream& out, const Edge& e) const {
        out << "[" << name[e].weight << "]";
        // out << "[weight=\"" << name[e].weight << "\"]";
    }
private:
    Name name;
};

template <class Name>
class custom_vertex_label_writer {
public:
    custom_vertex_label_writer(Name _name) : name(_name) {}
    template <class Vertex>
    void operator()(std::ostream& out, const Vertex& v) const {
        auto seq = name[v].read | seqan3::views::to_char;
        out << " [" << name[v].count << ", " << std::string(seq.begin(), seq.end()) << "]";
    }
private:
    Name name;
};

// Custom DFS visitor with edge condition and vertex saving
struct MyDFSVisitor : default_dfs_visitor {
    MyDFSVisitor(Graph& g, Vertex start_vertex, int max_L)
        : graph(g), startVertex(start_vertex), remainingValue(max_L) {}

    template <typename Vertex, typename Graph>
    void discover_vertex(Vertex u, const Graph& g) const {
        if (u != startVertex) {
            // Check if remaining value allows visiting this vertex
            for (auto it = adjacent_vertices(u, g); it.first != it.second; ++it.first) {
                auto v = *it.first;
                Edge e;
                bool exists;
                tie(e, exists) = edge(u, v, g);
                if (exists) {
                    EdgeProperties edge_props = g[e];
                    int edge_weight = edge_props.weight;

                    // Check if visiting this neighbor does not exceed remaining value
                    if (remainingValue - edge_weight >= 0) {
                        // Save the vertex as a potential new edge
                        savedVertices.insert(v);
                    }
                }
            }
        }
    }

    std::unordered_set<Vertex> savedVertices; // Container to store saved vertices

private:
    Graph& graph;
    Vertex startVertex;
    mutable uint8_t remainingValue;
};


class GraphConstructor
{
public:
    // GraphConstructor(std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads, cmd_arguments args);
    GraphConstructor(std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads, std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    void init_graph();
    void insert_edge(std::vector<seqan3::dna5> read1, std::vector<seqan3::dna5> read2, int edit_dis);
    std::vector<std::vector<seqan3::dna5>> mergeUniqueReads(const std::vector<std::vector<std::vector<seqan3::dna5>>>& read_vectors);
    void construct_graph();
    std::unordered_set<Vertex> visitNeighbors(const Graph& g, Vertex start_vertex, uint8_t dis);
    void save_graph() const;
    void edge_summary();
    // void process_block();
    // void process_block(const std::vector<std::vector<seqan3::dna5>> &reads_vec, int min_s, int max_s);
    // std::map<std::set<std::vector<seqan3::dna5>>, int> minimizer_omh();
    // std::map<std::set<std::vector<seqan3::dna5>>, int> omh_minimizer();
    // std::map<std::set<std::vector<seqan3::dna5>>, int> edges_main();
    // std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> get_edge_lst();
    // std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> get_edge_lst();
    // void display_edge_summary(std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst);
    // std::vector<std::pair<uint64_t, uint64_t>> get_combinations(const std::vector<uint64_t>& a);
    // std::unordered_set<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, customHash> unique_combination();
private:
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads_;
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2reads_;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads_;
    std::map<std::vector<seqan3::dna5>, uint32_t> read2count_;
    cmd_arguments args;
    std::filesystem::path graph_full_path_;
    std::unordered_map<std::vector<seqan3::dna5>, boost::graph_traits<Graph>::vertex_descriptor, std::hash<std::vector<seqan3::dna5>>> read2vertex_;
    std::unordered_map<boost::graph_traits<Graph>::vertex_descriptor, std::vector<seqan3::dna5>, std::hash<boost::graph_traits<Graph>::vertex_descriptor>> vertex2read_;
    Graph graph_;

    // std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> edge_lst;
    // std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> edge_lst;
    // std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst;    
    // std::map<int, size_t> edit_distance_counts_;
};
#endif // GraphConstructor_H



// Define a custom comparator for pairs that considers the order of elements
// struct pair_comparator
// {
//     template <typename T1, typename T2>
//     bool operator()(const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs) const
//     {
//         return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
//     }
// };

// Improved unordered_pair using boost::hash_combine
// struct customHash {
//     template <class T1, class T2>
//     std::size_t operator () (const std::pair<T1, T2> &p) const {
//         std::size_t seed = 0;
//         auto seq1 = p.first | seqan3::views::to_char;
//         string seq_str1(seq1.begin(), seq1.end());
//         auto seq2 = p.second | seqan3::views::to_char;
//         string seq_str2(seq2.begin(), seq2.end());
//         // compare 
//         // // auto seq1_hash = boost::hash_combine(static_cast<std::size_t>(1), seq_str1);
//         // // auto seq2_hash = boost::hash_combine(static_cast<std::size_t>(1), seq_str2);
//         // std::size_t seq1_hash = std::hash<std::string>{}(seq_str1);
//         // std::size_t seq2_hash = std::hash<std::string>{}(seq_str2);
//         boost::hash_combine(seed, seq_str1);
//         boost::hash_combine(seed, seq_str2);
//         return seed;
//     }
// };