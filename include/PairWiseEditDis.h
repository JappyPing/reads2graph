// PairWiseEditDis.h
#ifndef PAIRWISEEDITDIS_H
#define PAIRWISEEDITDIS_H

#include "Utils.h"
#include "LoggingLevels.h"
#include "EdgeConstructor.h"
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
    PairWiseEditDis(std::set<std::vector<seqan3::dna5>> const & reads, cmd_arguments args);

    std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> compute_pairwise_edit_distance();

private:
    std::set<std::vector<seqan3::dna5>> const & reads_;
    cmd_arguments args;
    std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> edge_lst;
    std::vector<std::vector<int>> pairwise_distances_;
    std::map<int, size_t> edit_distance_counts_;
};
#endif //PAIRWISEEDITDIS_H