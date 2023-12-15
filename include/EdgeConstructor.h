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
#include "Utils.h"
#include "LoggingLevels.h"

using namespace std;
using namespace seqan3::literals;

// Define a custom comparator for pairs that considers the order of elements
struct pair_comparator
{
    template <typename T1, typename T2>
    bool operator()(const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs) const
    {
        return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
    }
};

class EdgeConstructor
{
public:
    EdgeConstructor(std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads, cmd_arguments args);
    // void process_block();
    void process_block(const std::vector<std::vector<seqan3::dna5>> &reads_vec, int min_s, int max_s);
    void process_blocks_in_parallel();
    std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> get_edge_lst();
private:
    std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads_;
    cmd_arguments args;
    std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> edge_lst;
    std::map<int, size_t> edit_distance_counts_;
};
#endif // EDGECONSTRUCTOR_H