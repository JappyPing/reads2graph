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
#include "GraphConstructor.hpp"
#include "OMH.hpp"

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

int main(int argc, char** argv) {
    // check log whether log file exist or not, if exist, remove it
    // char cwd[1024];
    // if (getcwd(cwd, sizeof(cwd)) != NULL) {
    //     std::string filePath = std::string(cwd) + "/reads2graph" + oss.str();
    //     if (std::filesystem::exists(filePath)) {
    //         std::filesystem::remove(filePath);
    //     }
    // } else {
    //     std::cerr << "Error: unable to get current working directory." << std::endl;
    // }
    // say hello
    // 
	// Utils::getInstance().logger(LOG_LEVEL_INFO,  "Welcome to use reads2graph!");
    Utils::getInstance().logger(LOG_LEVEL_INFO,  "Welcome to use reads2graph!");

    std::string concatenatedArgs;

    for (int i = 0; i < argc; ++i) {
        concatenatedArgs += argv[i];
        if (i < argc - 1) {
            concatenatedArgs += " ";
        }
    }
    sharg::parser reads2graphParser{"reads2graph", argc, argv, sharg::update_notifications::off}; // initialise parser
    cmd_arguments args{};

    // utils.initialise_parser(reads2graphParser, args);
    Utils::getInstance().initialise_parser(reads2graphParser, args);

    Utils::getInstance().logger(LOG_LEVEL_INFO, concatenatedArgs);

    try
    {
        reads2graphParser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << "[Invalid Options] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Parameters: -o {} -k {} -w {} --omh_kmer_n {} --omh_times {}", args.output_dir.string(), args.k_size, args.window_size, args.omh_kmer_n, args.omh_times));

    // Declare and define a global variable for available cores
    int available_cores = omp_get_max_threads();
    Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("The maximum number of threads available: {} ", available_cores));
    // Ensure the user-specified number of cores is within a valid range
    auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    // std::cout << "The maximum number of threads available:" << available_cores << std::endl;
    
    // Set the number of threads for OpenMP
    // omp_set_num_threads(num_cores_to_use);
    args.num_process = num_cores_to_use;
    // std::cout << "The number of threads :" << num_cores_to_use << std::endl;
    Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of threads: {} ", num_cores_to_use));
    ////////////////////////////////////////////////////////////////////////////
    // auto read2count = ReadWrite(args).get_unique_reads_counts();
    auto [unique_reads, read2count, min_read_length] = ReadWrite(args).get_unique_reads_counts();
    args.read_length = min_read_length;
    auto total_uniq_num = unique_reads.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of unique reads: {}, minimum read length: {}.", total_uniq_num, min_read_length));

    // auto uniq_num = unique_reads.size();
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of unique reads: {} ", uniq_num));

    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    // #pragma omp parallel
    // {
    //     std::vector<std::vector<seqan3::dna5>> private_unique_reads;

    //     #pragma omp for nowait
    //     for (long unsigned int i = 0; i < uniq_num; ++i) {
    //         auto it = read2count.begin();
    //         std::advance(it, i);
    //         private_unique_reads.push_back(it->first);
    //     }

    //     #pragma omp critical
    //     unique_reads.insert(unique_reads.end(), private_unique_reads.begin(), private_unique_reads.end());
    // }

    // auto [read2count, read2id] = ReadWrite(args).get_unique_reads_counts();
    // // Print the unique reads.
    // for (auto const & read : unique_reads)
    // {
    //     seqan3::debug_stream << read << '\n';
    // }
    std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst;
    if (args.pair_wise) {
        GraphConstructor graph_constructor(read2count, args);
        graph_constructor.construt_graph_via_pairwise_comparison(unique_reads);
        if (args.save_graph){
            graph_constructor.save_graph();
        }
        // Create an instance of PairwiseEditDistance
        // edge_lst = PairWiseEditDis(unique_reads, args).compute_pairwise_edit_distance();
        // Utils::getInstance().logger(LOG_LEVEL_INFO,  "Pairwise (brute force) edit-distance-based edges calculation done");
        // GraphManager(edge_lst, read2count, args).construct_graph();
        // GraphManager(edge_lst, read2count, read2id, args).construct_graph();
    } else {
        // minimizer grouping first and then omh
        auto betterParams = MinimizerGenerator(args).possibleBetterParameters();

        auto hash2reads = MinimizerGenerator(args).minimizer2reads_main(unique_reads, betterParams);  

        GraphConstructor graph_constructor(read2count, args);
        graph_constructor.construct_graph(hash2reads);
        if (args.save_graph){
            graph_constructor.save_graph();
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  "edit-distance-based graph construction done!");
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