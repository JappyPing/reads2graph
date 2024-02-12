#ifndef READWRITE_HPP
#define READWRITE_HPP

#include "Utils.hpp"
#include <array>  // std::array
#include <string> // std::string
#include <vector> // std::vector
 
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/utility/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/io/sequence_file/all.hpp>

using namespace std;

class ReadWrite{
    public:
        ReadWrite(cmd_arguments args);//constructor
        ~ReadWrite(); //deconstructor
        // std::pair<std::set<std::vector<seqan3::dna5>>, std::map<std::vector<seqan3::dna5>, uint32_t>> get_unique_reads_counts(cmd_arguments args);
        std::map<std::vector<seqan3::dna5>, uint32_t> get_unique_reads_counts();
        // std::tuple<std::map<std::vector<seqan3::dna5>, uint32_t>, std::unordered_map<std::string, std::vector<seqan3::dna5>> get_unique_reads_counts();
    private:
        cmd_arguments args;
};

#endif