// OMH.h
#ifndef __OMH_H__
#define __OMH_H__

#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <limits>

#include <xxhash.hpp>
#include <seeded_prg.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

struct mer_info {
  size_t pos;
  uint64_t hash;
  unsigned occ;
  mer_info(size_t p, unsigned o, uint64_t h)
    : pos(p)
    , hash(h)
    , occ(o)
  { }
};

// template<typename EngineT>
// std::vector<uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned m, EngineT& prg);

template<typename EngineT>
// std::vector<uint64_t> omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, EngineT& prg) {
std::vector<uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned m, EngineT& prg) {
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
    std::vector<uint64_t> omh_results;
    for(unsigned i = 0; i < m; ++i) {
        const auto seed = prg();
        for(auto& meri : mers) {
        hash.reset(seed);
        hash.update(&seq.data()[meri.pos], k);
        if(weight) hash.update(&meri.occ, sizeof(meri.occ));
        meri.hash = hash.digest();
        }

        // std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), & { return x.hash < y.hash; });
        // std::partial_sort(mers.begin(), mers.begin() + l, mers.end(),  { return x.hash < y.hash; });

        std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
        std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });


        omh_results.push_back(mers[0].hash); // Save the smallest hash as the OMH result for this round
    }

    return omh_results; // Return the OMH results
}


#endif /* __OMH_H__ */
