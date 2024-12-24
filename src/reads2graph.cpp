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
#include "gOMH.hpp"

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <seqan3/core/debug_stream.hpp> // for debug_stream
// #include <format>
#include <boost/format.hpp>

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
    bool num_process_input = false;
    for (int i = 0; i < argc; ++i) {
        concatenatedArgs += argv[i];
        if (i < argc - 1) {
            concatenatedArgs += " ";
        }
        if (strcmp(argv[i], "--num_process") == 0 || strcmp(argv[i], "-p") == 0) {
            num_process_input = true;
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
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Parameters: -o {} -k {} -w {} --gomh_kmer_n {} --gomh_times {}", args.output_dir.string(), args.k_size, args.window_size, args.gomh_kmer_n, args.gomh_times));
    omp_set_dynamic(0);
    setenv("OMP_PROC_BIND", "false", 1);
    omp_set_num_threads(args.num_process);
    // Determine the number of cores available
    int available_cores = omp_get_max_threads();
    // Log the maximum number of threads available
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The maximum number of CPU cores available: %1%") % available_cores));

    // Ensure the number of cores to use is within a valid range
    if (num_process_input){
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of CPU cores requested: %1% ") % args.num_process));
        if (args.num_process > available_cores){
            args.num_process = available_cores;
        }
    } 
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of CPU cores actually used: %1% ") % args.num_process));
    ///////////////////////////////////////////////////
    // Declare and define a global variable for available cores
    // int available_cores = omp_get_max_threads();
    // // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("The maximum number of threads available: {} ", available_cores));
    // Utils::getInstance().logger(LOG_LEVEL_DEBUG, boost::str(boost::format("The maximum number of threads available: %1% ") % available_cores));
    // // Ensure the user-specified number of cores is within a valid range
    // auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    // std::cout << "The maximum number of threads available:" << available_cores << std::endl;
    
    // Set the number of threads for OpenMP
    // omp_set_num_threads(num_cores_to_use);
    // args.num_process = num_cores_to_use;
    // std::cout << "The number of threads :" << num_cores_to_use << std::endl;
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of threads: {} ", num_cores_to_use));

    ////////////////////////////////////////////////////////////////////////////
    // auto read2count = ReadWrite(args).get_unique_reads_counts();
    auto [unique_reads, read2count, min_read_length] = ReadWrite(args).get_unique_reads_counts();
    args.read_length = min_read_length;
    auto total_uniq_num = unique_reads.size();

    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of unique reads: %1%, minimum read length: %2%.") % total_uniq_num % min_read_length));

    // Validate the input mode
    if (!is_valid_bucketing_mode(args.bucketing_mode)) {
        std::cerr << "Error: Invalid bucketing mode selected! Choose from: minimizer_gomh, miniception, omh, brute_force. Default minimizer_gomh. miniception, omh and brute_force are implemented to assess the performance of reads2graph." << std::endl;
        return EXIT_FAILURE;
    }

    std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst;
    GraphConstructor graph_constructor(read2count, args);
    if (args.bucketing_mode == "brute_force") {
        graph_constructor.construt_graph_via_pairwise_comparison(unique_reads);
    } else if (args.bucketing_mode == "omh") {
        graph_constructor.construt_graph_via_omh(unique_reads);
    } else if (args.bucketing_mode == "miniception") {
        graph_constructor.construt_graph_via_miniception(unique_reads);   
    } else if (args.bucketing_mode == "minimizer_gomh") {
        auto hash2reads = MinimizerGenerator(args).minimizer2reads_main(unique_reads);
        graph_constructor.construct_graph(hash2reads);
    }
    if (args.save_graph){
        graph_constructor.save_graph();
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