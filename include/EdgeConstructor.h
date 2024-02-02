// EdgeConstructor.h
#ifndef EDGECONSTRUCTOR_H
#define EDGECONSTRUCTOR_H

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/all.hpp>
#include <vector>
#include <set>
#include <seqan3/alphabet/all.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits> // for std::decay_t
#include <boost/functional/hash.hpp>
#include "Utils.h"
#include "LoggingLevels.h"
#include <unordered_set>

using namespace std;
using namespace seqan3::literals;

// Define a custom comparator for pairs that considers the order of elements
// struct pair_comparator
// {
//     template <typename T1, typename T2>
//     bool operator()(const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs) const
//     {
//         return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
//     }
// };

// Improved unordered_pair using boost::hash_combine
// struct unordered_pair {
//     template <class T1, class T2>
//     std::size_t operator () (const std::pair<T1, T2> &p) const {
//         std::size_t seed = 0;
//         auto seq1 = p.first | seqan3::views::to_char;
//         string seq_str1(seq1.begin(), seq1.end());
//         auto seq2 = p.second | seqan3::views::to_char;
//         string seq_str2(seq2.begin(), seq2.end());
//         // compare 
//         // auto seq1_hash = boost::hash_combine(static_cast<std::size_t>(1), seq_str1);
//         // auto seq2_hash = boost::hash_combine(static_cast<std::size_t>(1), seq_str2);
//         std::size_t seq1_hash = std::hash<std::string>{}(seq_str1);
//         std::size_t seq2_hash = std::hash<std::string>{}(seq_str2);
//         if (seq1_hash == seq2_hash) {
//             // throw std::runtime_error("Error: two reads are the same in the unoredered read pair!");
//             std::cout << "The same!" << endl;
//         } else if (seq1_hash < seq2_hash) {
//             boost::hash_combine(seed, seq_str1);
//             boost::hash_combine(seed, seq_str2);
//         } else {
//             boost::hash_combine(seed, seq_str2);
//             boost::hash_combine(seed, seq_str1);        
//         }
//         return seed;
//     }
// };

class EdgeConstructor
{
public:
    EdgeConstructor(std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads, cmd_arguments args);
    // void process_block();
    void process_block(const std::vector<std::vector<seqan3::dna5>> &reads_vec, int min_s, int max_s);
    void process_blocks_in_parallel();
    // std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> get_edge_lst();
    // std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> get_edge_lst();
    std::map<std::set<std::vector<seqan3::dna5>>, int> get_edge_lst();
private:
    std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads_;
    cmd_arguments args;
    // std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> edge_lst;
    // std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> edge_lst;
    std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst;    
    std::map<int, size_t> edit_distance_counts_;
};
#endif // EDGECONSTRUCTOR_H