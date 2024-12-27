// MinimizerGenerator.hpp

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <vector>
#include <set>
#include <seqan3/alphabet/all.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits> // for std::decay_t
#include "Utils.hpp"
#include "LoggingLevels.hpp"

using namespace std;
using namespace seqan3::literals;

class MinimizerGenerator
{
public:
    MinimizerGenerator(cmd_arguments args);
    // MinimizerGenerator(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::tuple<unsigned, unsigned, unsigned, double> betterParams);
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads);
    double proba(unsigned L, unsigned k);
    int kSize(int L, double p);
    std::vector<std::vector<seqan3::dna5>> divide_into_substrings(const std::vector<seqan3::dna5> & read, int num_substrs);
    uint8_t k_estimate(uint8_t num_substr, uint8_t substr_size);
    uint8_t wSize(uint8_t k);
    // std::tuple<unsigned, unsigned, unsigned, double> overlappingWindowParameters();
    // std::tuple<unsigned, unsigned, unsigned, double> segmentation_parameters();
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer_only2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, uint8_t k, unsigned m);
    // std::tuple<int, int, double> findBestParameters(int l, int dt, double pt);
    // long double prob(int l, int n, int k, int dt);
    // void process_read(const std::vector<seqan3::dna5> &read);
    // std::size_t min_k_in_sampling(const std::vector<std::vector<seqan3::dna5>>& unique_reads, double p);

private:
    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads;
    cmd_arguments args;
};