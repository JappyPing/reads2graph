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

// Define a custom function to remove edges with weights in the specified interval
template <typename Graph>
void remove_edges_in_interval(Graph& g, int lower_bound, int upper_bound) {
    // Iterate over all edges in the graph
    typename graph_traits<Graph>::edge_iterator ei, ei_end, next;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei = next) {
        next = ei;
        ++next;
        // Check if the weight of the current edge falls within the specified interval
        if (g[*ei].weight >= lower_bound && g[*ei].weight <= upper_bound) {
            // If the weight is within the interval, remove the edge
            remove_edge(*ei, g);
        }
    }
}

class GraphConstructor
{
public:
    GraphConstructor(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    void init_graph();
    void insert_edge(std::vector<seqan3::dna5> read1, std::vector<seqan3::dna5> read2, int edit_dis);
    std::vector<std::vector<seqan3::dna5>> mergeUniqueReads(const std::vector<std::vector<std::vector<seqan3::dna5>>>& read_vectors);
    void construct_graph(std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads);

    void visitNeighborsWithThreshold(const Graph& g, Vertex node, int distance_threshold, int current_distance, std::vector<Vertex>& indirect_neighbors, std::vector<bool>& visited);
    std::vector<Vertex> visitNeighborsOfNeighborsWithThreshold(const Graph& g, Vertex node, int distance_threshold);
    void save_graph() const;
    void edge_summary();
    void construt_graph_via_pairwise_comparison(std::vector<std::vector<seqan3::dna5>> unique_reads);

private:
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads_;
    std::map<std::vector<seqan3::dna5>, uint32_t> read2count_;
    cmd_arguments args;
    std::filesystem::path graph_full_path_;
    std::unordered_map<std::vector<seqan3::dna5>, Vertex, std::hash<std::vector<seqan3::dna5>>> read2vertex_;
    std::unordered_map<Vertex, std::vector<seqan3::dna5>, std::hash<Vertex>> vertex2read_;
    Graph graph_;
};
#endif // GraphConstructor_H
