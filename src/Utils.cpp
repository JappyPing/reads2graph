/**
 * @Author: Pengyao Ping
 * @Date:   2023-02-07 12:51:41
 * @Last Modified by:   Pengyao Ping
 * @Last Modified time: 2023-03-05 12:47:18
 */

#include "Utils.h"
#include "LoggingLevels.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include <unordered_map>
#include <set>
#include <map>
#include <thread>
#include <unistd.h>
#include <functional>
// #include <pbbam/DataSet.h>


#define ReadGraph_VERSION "0.1.1"

using namespace std;
// using namespace PacBio;
// using namespace BAM;

Utils::Utils(void){

}

Utils::~Utils(void){

}

// Define a logging function to log messages to a log file when necessary
void Utils::logMessage(int log_level, const std::string& message) {
    if (log_level < LOG_LEVEL) return;

    // Get the current time and format it as a string
    time_t rawTime;
    time(&rawTime);
    struct tm * timeInfo = localtime(&rawTime);
    char timeString[20];
    strftime(timeString, 20, "%Y-%m-%d %H:%M:%S", timeInfo);

    char cwd[1024];
    std::ofstream logFile;
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::string filePath = std::string(cwd) + "/ReadGraph.log";
        logFile.open(filePath, std::ios::app);
    } else {
        std::cerr << timeString << ": "<< "Error: unable to get current working directory." << std::endl;
    }
    
    std::string log_color;
    std::string log_prefix;
    switch (log_level) {
    case LOG_LEVEL_DEBUG:
        log_color = COLOR_DEBUG;
        log_prefix = "DEBUG: ";
        break;
    case LOG_LEVEL_INFO:
        log_color = COLOR_INFO;
        log_prefix = "INFO: ";
        break;
    case LOG_LEVEL_WARNING:
        log_color = COLOR_WARNING;
        log_prefix = "WARNING: ";
        break;
    case LOG_LEVEL_ERROR:
        log_color = COLOR_ERROR;
        log_prefix = "ERROR: ";
        break;
    case LOG_LEVEL_CRITICAL:
        log_color = COLOR_CRITICAL;
        log_prefix = "CRITICAL ERROR: ";
        break;
    default:
        log_color = COLOR_RESET;
        break;
    }

    std::cout << timeString << ": " << log_color << log_prefix << message << COLOR_RESET << std::endl;
    if(logFile.is_open()){
        logFile << timeString << ": " << log_prefix << message << endl;
        logFile.close();
    }else {
        std::cerr << "\033[1;31mERROR: ReadGraph log file opened failed.\033[0m" << std::endl;
    }
    
}
/////////////////////////////////////////////////////////////////

// Define a function to control all the parameters via command line
void Utils::initialise_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Pengyao Ping";
    parser.info.short_description = "Construct a read graph from a short read set.";
    parser.info.version = ReadGraph_VERSION;
 
    parser.add_option(args.input_data, 
                        sharg::config{.short_id = 'i', 
                                      .long_id = "input_data", 
                                      .description = "Please provide a fasta/fastq/ data file."});

    parser.add_option(args.output_dir,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output_dir",
                                    .description = "The directory for outputs."});

    parser.add_option(args.k_size,
                      sharg::config{.short_id = 'k',
                                    .long_id = "k_size",
                                    .description = "The size for minimiser."});

    parser.add_option(args.window_size,
                      sharg::config{.short_id = 'w',
                                    .long_id = "window_size",
                                    .description = "The window size for minimiser."});


    parser.add_option(args.max_edit_dis,
                      sharg::config{.short_id = 'x',
                                    .long_id = "max_edit_dis",
                                    .description = "The maximum edit distance for constructing edges between reads"});

    parser.add_option(args.min_edit_dis,
                      sharg::config{.short_id = 'n',
                                    .long_id = "min_edit_dis",
                                    .description = "The minimum edit distance for constructing edges between reads."});

    parser.add_option(args.num_process,
                      sharg::config{.short_id = 'p',
                                    .long_id = "num_process",
                                    .description = "The number of expected processes."});

    parser.add_option(args.graph_filename,
                      sharg::config{.short_id = 'g',
                                    .long_id = "graph_filename",
                                    .description = "The file name of the constructed graph."});
    // Brute Force
    parser.add_option(args.pair_wise,
                      sharg::config{
                                    .long_id = "pair_wise",
                                    .description = "Brute Force calcualte the pairwise edit distance."});

    // parser.add_option(args.year,
    //                   sharg::config{.short_id = 'y',
    //                                 .long_id = "year",
    //                                 .description = "Only data entries that are newer than `year` are considered."});
 
    // parser.add_option(args.aggregate_by,
    //                   sharg::config{.short_id = 'a',
    //                                 .long_id = "aggregate-by",
    //                                 .description = "Choose your method of aggregation: mean or median."});
 
    // parser.add_flag(
    //     args.header_is_set,
    //     sharg::config{.short_id = 'H',
    //                   .long_id = "header-is-set",
    //                   .description =
    //                       "Let us know whether your data file contains a header to ensure correct parsing."});
}




// // Define a function to control all the parameters via command line
// std::map<std::string, std::string> Utils::parseArgvs(int argc, char* argv[]) {
//     if (argc <=1 ) {
//         logMessage(LOG_LEVEL_ERROR, "No parameters were input.");
//         exit(1);
//     }
//     else{
//         std::map<std::string, std::string> paras;
//         for (int i = 1; i < argc; i++) {
//             std::string argument = argv[i];
//             if (argument == "-h" || argument == "--help") {
//                 printHelpMessage();
//                 exit(0);
//             } else if (argument == "-v" || argument == "--version"){
//                 std::cout << "ReadGraph version: " << ReadGraph_VERSION << std::endl;
//             }
//             else{
//                 std::string key;
//                 if (argument[0] == '-'){
//                     key = argument.substr(1);
//                 }else if (argument.substr(0, 2) == "--"){
//                     key = argument.substr(2);
//                 }

//                 if (stringExists(key)){
//                     if (i + 1 < argc && argv[i + 1][0] != '-') {
//                         paras[key] = argv[i + 1];
//                         i++;
//                     } 
//                     else {
//                         paras[key] = "";
//                     }    
//                 }
//                 else {
//                     logMessage(LOG_LEVEL_ERROR, "Invalid parameters were input");
//                     exit(1);
//                 }
//             }
//         }
//         return paras;
//     }
//     return std::map<std::string, std::string>();

// }

// bool Utils::stringExists(const std::string& input) {
//     std::set<std::string> parasAllowedValues = {"h", "help", "i", "input", "o", "output", "t", "thread"};
//     return parasAllowedValues.count(input) > 0;
// }

// // print help message
// void Utils::printHelpMessage() {
//     std::cout << "Usage: ReadGraph [OPTIONS]\n"
//             << "\n"
//             << "Options:\n"
//             << "  -h, --help            Show this help message and exit\n"
//             << "  -i, --input filename  Input filename\n"
//             << "  -o, --output filename Output filename\n"
//             << "  -v, --version output ReadGraph version\n"
//             << "  -t, --thread output ReadGraph version\n"
//             << std::endl;
// }

// // Get the number of hardware concurrency levels, i.e. the number of CPU cores
// int Utils::getHardwareConcurrency() {
//     return std::thread::hardware_concurrency();
// }

// std::vector<int> Utils::get_sub_lst(const std::vector<int>& input, int startIndex, int endIndex) {
//     return std::vector<int>(input.begin() + startIndex, input.begin() + endIndex);
// }

// std::vector<std::vector<int>> Utils::split_vector(const std::vector<int>& vec, int num_groups) {
//   std::vector<std::vector<int>> groups;
//   groups.resize(num_groups);

//   int group_size = vec.size() / num_groups;
//   int remaining = vec.size() % num_groups;

//   for (int i = 0, start = 0; i < num_groups; ++i) {
//     int size = group_size + (remaining > 0 ? 1 : 0);
//     if (remaining > 0) --remaining;

//     groups[i].assign(vec.begin() + start, vec.begin() + start + size);
//     start += size;
//   }

//   return groups;
// }




