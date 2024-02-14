#include "GraphManager.hpp"
#include "ReadWrite.hpp"
#include <iostream>

// GraphManager::GraphManager(std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> edge_lst, std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : edge_lst_(std::move(edge_lst)), read2count_(std::move(read2count)), args(args) {}

// GraphManager::GraphManager(std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> edge_lst, std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : edge_lst_(std::move(edge_lst)), read2count(read2count), args(args) {}

GraphManager::GraphManager(std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst, std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : edge_lst_(std::move(edge_lst)), read2count(read2count), args(args) {}
GraphManager::~GraphManager(){}

void GraphManager::construct_graph()
{
    for (const auto& [read, count] : read2count) {
        auto v = boost::add_vertex({read, count}, graph);
        read2vertex[read] = v;
        vertex2read[v] = read;
    }
    for (const auto& [read_pair_set, edit_distance] : edge_lst_) {
        auto it = read_pair_set.begin(); 
        std::vector<seqan3::dna5> first_read = *it;       
        auto v1_iter = read2vertex[first_read];
        ++it;
        std::vector<seqan3::dna5> second_read = *it;
        auto v2_iter = read2vertex[second_read];
        // auto v1_iter = read2vertex[read_pair.first];
        // auto v2_iter = read2vertex[read_pair.second];
        boost::add_edge(v1_iter, v2_iter, {first_read, second_read, edit_distance}, graph);
    }
    save_graph();
    // auto graph_full_path_ = args.output_dir / args.graph_filename;
    // Graph originalGraph = read_graph_from_dot(graph_full_path_);
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
    // auto graph_full_path_ = args.output_dir / args.graph_filename;
    std::filesystem::path modifiable_path = args.input_data;
    modifiable_path.replace_extension(".dot");
    std::string output_filename = modifiable_path.filename().string();
    auto graph_full_path_ = args.output_dir / output_filename;
    std::cout << "Output graph Path: " << graph_full_path_ << std::endl;
    // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Output graph Path: {}.", graph_full_path_));
    std::ofstream dot_file(graph_full_path_);

    // boost::write_graphviz(dot_file, graph);  
    
    // Check if the file is open
    if (dot_file.is_open()) {
        // Write the graph to the file in Graphviz format
        // boost::write_graphviz(dot_file, graph, make_label_writer(get(&VertexProperties::count, graph)),
        //                make_label_writer(get(&EdgeProperties::weight, graph)));

        // boost::write_graphviz(dot_file, graph, make_label_writer(get(&VertexProperties::count, graph)), custom_edge_label_writer<Graph>(graph));
        boost::write_graphviz(dot_file, graph, custom_vertex_label_writer<Graph>(graph), custom_edge_label_writer<Graph>(graph));

        // Close the file
        dot_file.close();
        
        std::cout << "Graph written to " << graph_full_path_ << " successfully." << std::endl;
        // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Graph written to {} successfully.", graph_full_path_));
    } else {
        std::cerr << "Error opening file " << graph_full_path_ << std::endl;
    }        
    // std::cout << "Graph saved to file successfully." << std::endl;
}

// Function to read the graph from a .dot file and create the original graph
// Graph GraphManager::read_graph_from_dot(const std::string& filename) {
//     std::ifstream file(filename);

//     if (!file.is_open()) {
//         std::cerr << "Error opening file " << filename << std::endl;
//         return Graph(); // Return an empty graph on error
//     }

//     Graph graph;

//     dynamic_properties dp;
//     dp.property("weight", get(&EdgeProperties::weight, graph));
//     dp.property("count", get(&VertexProperties::count, graph));
//     dp.property("read", get(&VertexProperties::read, graph));

//     read_graphviz(file, graph, dp);

//     file.close();
//     return graph;
// }