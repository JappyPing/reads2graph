/**
 * @Author: Pengyao Ping
 * @Date:   2023-02-07 14:10:30
 * @Last Modified by:   Pengyao Ping
 * @Last Modified time: 2023-03-05 12:52:34
 */
#ifndef UTILS_HPP
#define UTILS_HPP

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

struct cmd_arguments
{
    std::filesystem::path input_data{};
    std::filesystem::path output_dir{std::filesystem::current_path()};
    unsigned read_length{};
    std::size_t total_uniq_num{};
    bool default_params{true};
    uint8_t k_size{};
    uint8_t w_size{};
    uint8_t substr_number{3};
    uint8_t max_edit_dis{2};
    uint8_t min_edit_dis{1};
    int num_process{1};
    unsigned int bin_size_max{10000};
    unsigned gomh_k{4};
    unsigned gomh_times{3};
    bool gomh_flag{false};
    std::uint64_t seed{2025};
    double differ_kmer_ratio{0.3};
    double probability{0.86};
    unsigned visit_depth{15};
    bool save_graph{false};
    std::string bucketing_mode{"minimizer_gomh"};
    bool segmentation{true};
    bool miniception_gomh{true};

    unsigned omh_l{2};
    unsigned omh_m{3};
    unsigned omh_k{18};
    unsigned miniception_w{19};
    unsigned miniception_k{18};    
};

// Utility function to validate the bucketing mode
bool is_valid_bucketing_mode(const std::string &mode);

class Utils{
    public:
        static Utils& getInstance();

        void logger(int log_level, const std::string& message);
        void initialise_parser(sharg::parser & parser, cmd_arguments & args);

        // Utils(void);//constructor
        // ~Utils(void); //deconstructor
        // static void logMessage(int log_level, const std::string& message);   
        // void initialise_parser(sharg::parser & parser, cmd_arguments & args);
        ////////////////////////////////////////////
        // static bool stringExists(const std::string& input);
        // static std::map<std::string, std::string> parseArgvs(int argc, char* argv[]);
        // static void printHelpMessage();
        // static int getHardwareConcurrency();
        // static std::vector<int> get_sub_lst(const std::vector<int>& input, int startIndex, int endIndex);
        // static std::vector<std::vector<int>> split_vector(const std::vector<int>& vec, int num_groups);
    private:
        std::ofstream logFile;

        Utils();
        Utils(const Utils&) = delete;
        Utils& operator=(const Utils&) = delete;        

};

#endif