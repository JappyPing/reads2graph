/*
 * @Author: Pengyao PING
 * @Date: 2021-07-19 17:24:08
 * @LastEditors: Pengyao PING
 * @LastEditTime: 2021-08-11 11:36:10
 * @Email: Pengyao.Ping@student.uts.edu.au
 * @Description: 
 */

#include "Utils.h"
#include "LoggingLevels.h"
#include "ReadWrite.h"
#include "MinimizerGenerator.h"
#include "EdgeConstructor.h"
#include "GraphManager.h"
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
#include <iostream>
#include <vector>

using namespace std;
using namespace seqan3::literals;
// global default parameters
// bool logToFile = true;
// int numThreads = 1;

int main(int argc, char** argv) {
    // check log whether log file exist or not, if exist, remove it
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::string filePath = std::string(cwd) + "/ReadGraph.log";
        if (std::filesystem::exists(filePath)) {
            std::filesystem::remove(filePath);
        }
    } else {
        std::cerr << "Error: unable to get current working directory." << std::endl;
    }
    // say hello
    // 
	Utils::logMessage(LOG_LEVEL_INFO,  "Welcome to use ReadGraph!");
    // seqan3::debug_stream << "Hello World!\n";

    sharg::parser ReadGraphParser{"ReadGraph", argc, argv}; // initialise parser
    cmd_arguments args{};
 
    Utils utils;
    utils.initialise_parser(ReadGraphParser, args);
 
    try
    {
        ReadGraphParser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << "[Invalid Options] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    ReadWrite readWrite;
    auto results = readWrite.get_unique_reads_counts(args);
    // // Print the unique reads.
    // for (auto const & read : unique_reads)
    // {
    //     seqan3::debug_stream << read << '\n';
    // }
    auto minimiser_to_reads = MinimizerGenerator(results.first, args).process_reads_in_parallel();
    // EdgeConstructor(minimiser_to_reads, args).process_block();   

    auto edge_lst = EdgeConstructor(minimiser_to_reads, args).get_edge_lst();
    ///////////////////////////////////////////////////////////////////////////////
    GraphManager(edge_lst, results.second, args).construct_graph();
    
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

    

    // const int k = 5;
    // std::string inputFilePath = "../data/demo.fa"; // Replace with your actual file path
    // std::string query = "ATCGATCAG";
    // std::string outputFilePath = "output_graph.txt"; // Replace with your desired output file path

    // WorkflowController workflowController;
    // workflowController.runWorkflow(inputFilePath, query, paras["k"], outputFilePath);


    // // split the large pbbam file into few files up to the number of cpu cores - 2
    // DataProcess DP;
    // int cpu_core_num = Utils::getHardwareConcurrency();
    // // DP.splitPBBamFile(paras["i"], cpu_core_num-2, outputFolder);
    // Controller main_controller;
    // main_controller.pbbam2records(paras["i"], std::stoi(paras["t"]));

    return 0;
}