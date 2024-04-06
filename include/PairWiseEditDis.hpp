// PairWiseEditDis.h
#ifndef PAIRWISEEDITDIS_HPP
#define PAIRWISEEDITDIS_HPP

#include "Utils.hpp"
#include "LoggingLevels.hpp"
// #include "EdgeConstructor.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/pairwise/all.hpp>

using namespace std;
using namespace seqan3::literals;

class PairWiseEditDis
{
public:
    // Constructor
    // PairWiseEditDis(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    PairWiseEditDis(std::vector<std::vector<seqan3::dna5>> unique_reads, cmd_arguments args);
    // std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> compute_pairwise_edit_distance();
    std::map<std::set<std::vector<seqan3::dna5>>, int> compute_pairwise_edit_distance();

private:
    // std::set<std::vector<seqan3::dna5>> const & reads_;
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    std::vector<std::vector<seqan3::dna5>> unique_reads;
    cmd_arguments args;
    // std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> edge_lst;
    std::map<std::set<std::vector<seqan3::dna5>>, int> edge_lst;
    std::vector<std::vector<int>> pairwise_distances_;
    std::map<int, size_t> edit_distance_counts_;
};
#endif //PAIRWISEEDITDIS_H