// ReadGraph.cpp

#include "ReadGraph.hpp"

// namespace reads2graph{
    ReadGraph::ReadGraph(cmd_arguments args) : args(args) {}

    Graph ReadGraph::reads_graph(){
        int available_cores = omp_get_max_threads();
        // Ensure the user-specified number of cores is within a valid range
        auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
        args.num_process = num_cores_to_use;
        Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of threads: {} ", num_cores_to_use));

        auto [unique_reads, read2count, min_read_length] = ReadWrite(args).get_unique_reads_counts();
        args.read_length = min_read_length;
        auto total_uniq_num = unique_reads.size();
        Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of unique reads: {} ", total_uniq_num));

        auto betterParams = MinimizerGenerator(args).possibleBetterParameters();
        args.omh_k = std::get<2>(betterParams);

        auto hash2reads = MinimizerGenerator(args).minimizer2reads_main(unique_reads, betterParams); 

        auto edge_lst = EdgeConstructor(hash2reads, args).edges_main(); 
        Utils::getInstance().logger(LOG_LEVEL_INFO,  "Edit-distance-based edges calculation done!");

        Graph graph;

        if (args.save_graph){
            GraphManager graph_manager(edge_lst, read2count, args);
            graph = graph_manager.construct_graph();
            graph_manager.save_graph();                       
        } else {
            graph = GraphManager(edge_lst, read2count, args).construct_graph(); 
        }
        return graph;
    }
// }
