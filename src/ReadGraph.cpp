/*
 * @Author: Pengyao PING
 * @Date: 2021-07-19 17:24:08
 * @LastEditors: Pengyao PING
 * @LastEditTime: 2021-08-11 11:36:10
 * @Email: Pengyao.Ping@student.uts.edu.au
 * @Description: 
 */

#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include "ReadWrite.hpp"
#include "MinimizerGenerator.hpp"
#include "EdgeConstructor.hpp"
#include "GraphManager.hpp"
#include "PairWiseEditDis.hpp"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <seqan3/core/debug_stream.hpp> // for debug_stream

#include <array>  // std::array
#include <string> // std::string
#include <vector> // std::vector
 
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/utility/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/io/sequence_file/all.hpp>
// #include <seqan3/std/filesystem>

#include <omp.h>

using namespace std;
using namespace seqan3::literals;
// global default parameters
// bool logToFile = true;
// int numThreads = 1;

int main(int argc, char** argv) {
    // check log whether log file exist or not, if exist, remove it
    // char cwd[1024];
    // if (getcwd(cwd, sizeof(cwd)) != NULL) {
    //     std::string filePath = std::string(cwd) + "/ReadGraph" + oss.str();
    //     if (std::filesystem::exists(filePath)) {
    //         std::filesystem::remove(filePath);
    //     }
    // } else {
    //     std::cerr << "Error: unable to get current working directory." << std::endl;
    // }
    // say hello
    // 
	// Utils::getInstance().logger(LOG_LEVEL_INFO,  "Welcome to use ReadGraph!");
    Utils::getInstance().logger(LOG_LEVEL_INFO,  "Welcome to use ReadGraph!");

    sharg::parser ReadGraphParser{"ReadGraph", argc, argv}; // initialise parser
    cmd_arguments args{};
 
    // Utils utils;
    // utils.initialise_parser(ReadGraphParser, args);
    Utils::getInstance().initialise_parser(ReadGraphParser, args);
 
    try
    {
        ReadGraphParser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << "[Invalid Options] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    // Declare and define a global variable for available cores
    int available_cores = omp_get_max_threads();
    // Ensure the user-specified number of cores is within a valid range
    auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    std::cout << "The maximum number of threads available:" << available_cores << std::endl;
    // Set the number of threads for OpenMP
    omp_set_num_threads(num_cores_to_use);
    std::cout << "The number of threads setted:" << num_cores_to_use << std::endl;
    ////////////////////////////////////////////////////////////////////////////
    auto read2count = ReadWrite(args).get_unique_reads_counts();
    // auto [read2count, read2id] = ReadWrite(args).get_unique_reads_counts();
    // // Print the unique reads.
    // for (auto const & read : unique_reads)
    // {
    //     seqan3::debug_stream << read << '\n';
    // }

    if (args.pair_wise) {
            // Create an instance of PairwiseEditDistance
            auto edge_lst = PairWiseEditDis(read2count, args).compute_pairwise_edit_distance();
            Utils::getInstance().logger(LOG_LEVEL_INFO,  "Pairwise (brute force) nt-edit-distance-based edges calculation done");
            GraphManager(edge_lst, read2count, args).construct_graph();
            // GraphManager(edge_lst, read2count, read2id, args).construct_graph();
    } else {
        auto minimiser_to_reads = MinimizerGenerator(read2count, args).process_reads_in_parallel();
        // EdgeConstructor(minimiser_to_reads, args).process_block();   

        auto edge_lst = EdgeConstructor(minimiser_to_reads, args).get_edge_lst();
        Utils::getInstance().logger(LOG_LEVEL_INFO,  "nt-edit-distance-based edges calculation done!");
        ///////////////////////////////////////////////////////////////////////////////
        
        GraphManager(edge_lst, read2count, args).construct_graph();
        // GraphManager(edge_lst, read2count, read2id, args).construct_graph();

    }
    Utils::getInstance().logger(LOG_LEVEL_INFO,  "nt-edit-distance-based graph construction done!");
    //Print the stored read pairs and edit distances
    // for (const auto &[read_pair, edit_distance] : edge_lst)
    // {
    //     seqan3::debug_stream << "Read Pair: " << read_pair.first << " - " << read_pair.second
    //                         << ", Edit Distance: " << edit_distance << '\n';
    // }   

    // Parse command line arguments
    // std::map<std::string, std::string> paras;
    // paras = Utils::parseArgvs(argc, argv);
	
    // // ouput directory
    // std::string outputFolder;
    // if (paras.count("o") > 0) {
    //     outputFolder = paras["o"];
    // } else if ((paras.count("output_dir") > 0)) {
    //     outputFolder = paras["output_dir"];
    // } else {
    //     outputFolder = std::string(cwd) + "/result/";
    //     if (!std::filesystem::exists(outputFolder)) {
    //         std::filesystem::create_directory(outputFolder);
    //     }        
    // }

    



    return 0;
}