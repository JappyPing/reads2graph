#include "GraphManager.h"
#include <iostream>

GraphManager::GraphManager(std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> edge_lst, std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : edge_lst_(std::move(edge_lst)), read2count_(std::move(read2count)), args(args) {}

GraphManager::~GraphManager(){}


void GraphManager::construct_graph()
{
    for (const auto& [read, count] : read2count_) {
        auto v = boost::add_vertex({read, count}, graph);
        read2vertex[read] = v;
        vertex2read[v] = read;
    }
    for (const auto& [read_pair, edit_distance] : edge_lst_) {
        auto v1_iter = read2vertex[read_pair.first];
        auto v2_iter = read2vertex[read_pair.second];

        boost::add_edge(v1_iter, v2_iter, {read_pair.first, read_pair.second, edit_distance}, graph);

    }
    save_graph();
    // Iterate over vertices and print their attributes
    // boost::graph_traits<Graph>::vertex_iterator vi, vend;
    // for (boost::tie(vi, vend) = boost::vertices(graph); vi != vend; ++vi) {
    //     std::cout << "Node ID: " << *vi << ", Count: " << graph[*vi].count << ", Read: ";
    //     for (const auto &base : graph[*vi].read) {
    //         std::cout << base;
    //     }
    //     std::cout << std::endl;
    // }

    // // Iterate over edges and print their attributes
    // boost::graph_traits<Graph>::edge_iterator ei, eend;
    // for (boost::tie(ei, eend) = boost::edges(graph); ei != eend; ++ei) {
    //     std::cout << "Edge weight: " << graph[*ei].weight << std::endl;
    // }    
}

void GraphManager::save_graph() const {
    // std::ofstream out(graph_full_path_);
    auto graph_full_path_ = args.output_dir / args.graph_filename;
    std::cout << "Output graph Path: " << graph_full_path_ << std::endl;
    std::ofstream dot_file(graph_full_path_);
    boost::write_graphviz(dot_file, graph);  
    std::cout << "Graph saved to file successfully." << std::endl;
}