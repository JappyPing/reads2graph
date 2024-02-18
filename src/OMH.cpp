#include "OMH.hpp"
#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <format>
#include <seqan3/alphabet/all.hpp>

// OMH::OMH(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : read2count(read2count), args(args) {}
OMH::OMH(std::vector<std::vector<seqan3::dna5>> unique_reads, cmd_arguments args) : unique_reads(unique_reads), args(args) {}

std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH::omh2read_main(){
    int available_cores = omp_get_max_threads();
    auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    omp_set_num_threads(num_cores_to_use);

    std::default_random_engine prg;
    vector<unsigned int> seeds;
    for(unsigned i = 0; i < args.omh_times; ++i) {
        seeds.push_back(prg());
    }
    // auto uniq_num = read2count.size();
    // auto uniq_num = unique_reads.size();
    // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("The number of unique reads: {} ", uniq_num));
    // #pragma omp parallel for
    // for (size_t i = 0; i < uniq_num; ++i) {
    //     auto it = std::next(read2count.begin(), i);
    //     const auto& [read, count] = *it;
    #pragma omp parallel for
    for (auto const & read : unique_reads){

        for(auto seed : seeds){
            auto omh_value = omh_pos(read, 25, args.omh_kmer_n, seed);   //args.k_size
            #pragma omp critical
            {
                omh2reads[omh_value].push_back(read);
            }               
        } 
        // if (i % 10000 == 0) {
        //     // std::cout << "Iteration: " << i << '\n';
        //     Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("{} unique reads for OMH done!", i));
        // }
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
    return omh2reads;            
}

uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned int seed) {
    auto seql = read | seqan3::views::to_char;
    string seq(seql.begin(), seql.end());
    if(seq.size() < k) return {};

    const bool weight = l > 0;
    if(l == 0) l = 1;

    std::vector<mer_info> mers;
    std::unordered_map<std::string, unsigned> occurrences;

    //  Create list of k-mers with occurrence numbers
    for(size_t i = 0; i < seq.size() - k + 1; ++i) {
        auto occ = occurrences[seq.substr(i, k)]++;
        mers.emplace_back(i, occ, (uint64_t)0);
    }

    xxhash hash;
    // uint64_t omh_result;

    for(auto& meri : mers) {
      hash.reset(seed);
      hash.update(&seq.data()[meri.pos], k);
      if(weight) hash.update(&meri.occ, sizeof(meri.occ));
      meri.hash = hash.digest();
    }
    std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
    std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });

    return mers[0].hash; // Return the OMH results
}