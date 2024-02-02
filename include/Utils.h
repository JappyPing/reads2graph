/**
 * @Author: Pengyao Ping
 * @Date:   2023-02-07 14:10:30
 * @Last Modified by:   Pengyao Ping
 * @Last Modified time: 2023-03-05 12:52:34
 */
#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <memory>
#include <map>
#include <iostream>
#include <stdint.h>
#include <getopt.h>
#include <vector>
#include <functional>
// #include <pbbam/DataSet.h>
#include <sharg/all.hpp> // includes all necessary headers
#pragma once
extern int num_cores_to_use;  // Declare the global variable

using namespace std;
// using namespace PacBio;
// using namespace BAM;

struct cmd_arguments
{
    std::filesystem::path input_data{};
    std::filesystem::path output_dir{std::filesystem::current_path()};
    uint8_t k_size{4};
    uint8_t window_size{8};
    uint8_t max_edit_dis{5};
    uint8_t min_edit_dis{1};
    int num_process{26};
    std::string graph_filename{"graph.dot"};
    bool pair_wise{false};
};

class Utils{
    public:
        Utils(void);//constructor
        ~Utils(void); //deconstructor
        static void logMessage(int log_level, const std::string& message);   
        void initialise_parser(sharg::parser & parser, cmd_arguments & args);
        // static bool stringExists(const std::string& input);
        // static std::map<std::string, std::string> parseArgvs(int argc, char* argv[]);
        // static void printHelpMessage();
        // static int getHardwareConcurrency();
        // static std::vector<int> get_sub_lst(const std::vector<int>& input, int startIndex, int endIndex);
        // static std::vector<std::vector<int>> split_vector(const std::vector<int>& vec, int num_groups);

};

#endif