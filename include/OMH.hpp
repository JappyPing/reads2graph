// OMH.hpp

/*The implementation of Order Min Hash (OMH) is modified based on the original implementation of OMH(https://github.com/Kingsford-Group/omhismb2019). If you want to use the relevant source codes in your project, please remember to cite the original work listed below.
Guillaume Marçais, Dan DeBlasio, Prashant Pandey, Carl Kingsford, Locality-sensitive hashing for the edit distance, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i127–i135, https://doi.org/10.1093/bioinformatics/btz354
*/

#ifndef __OMH_HPP__
#define __OMH_HPP__

#include "Utils.hpp"
#include "LoggingLevels.hpp"

#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <limits>

// #include <xxhash.hpp>
// #include <seeded_prg.hpp>
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

class OMH
{
public:
    // OMH(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    OMH(cmd_arguments args);
    std::vector<std::pair<std::uint64_t, unsigned>> get_seeds_k();
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k);
    uint64_t omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed);
    
    unsigned omh_k(unsigned L, double p, uint8_t d);
    std::string getGappedSubstring(const std::string& str, size_t startPos, size_t length);
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::uint64_t seed, unsigned k);
    ////////////////////////////////
    // std::vector<uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed, int m);
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omhs2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned k, std::uint64_t seed, int m);
    ////////////////////////////////
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::uint64_t> seeds, unsigned k);
    
    // uint64_t omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned int seed);
private:
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2reads;
    cmd_arguments args;

};

// uint64_t omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned int seed) {
//     auto seql = read | seqan3::views::to_char;
//     string seq(seql.begin(), seql.end());
//     if(seq.size() < k) return {};

//     const bool weight = l > 0;
//     if(l == 0) l = 1;

//     std::vector<mer_info> mers;
//     std::unordered_map<std::string, unsigned> occurrences;

//     //  Create list of k-mers with occurrence numbers
//     for(size_t i = 0; i < seq.size() - k + 1; ++i) {
//         auto occ = occurrences[seq.substr(i, k)]++;
//         mers.emplace_back(i, occ, (uint64_t)0);
//     }

//     xxhash hash;
//     // uint64_t omh_result;

//     for(auto& meri : mers) {
//       hash.reset(seed);
//       hash.update(&seq.data()[meri.pos], k);
//       if(weight) hash.update(&meri.occ, sizeof(meri.occ));
//       meri.hash = hash.digest();
//     }
//     std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
//     std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });

//     return mers[0].hash; // Return the OMH results
// }


// template<typename EngineT>
// std::vector<uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned m, EngineT& prg);

// template<typename EngineT>
// // std::vector<uint64_t> omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, EngineT& prg) {
// std::vector<uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned m, EngineT& prg) {
//     auto seql = read | seqan3::views::to_char;
//     string seq(seql.begin(), seql.end());
//     if(seq.size() < k) return {};

//     const bool weight = l > 0;
//     if(l == 0) l = 1;

//     std::vector<mer_info> mers;
//     std::unordered_map<std::string, unsigned> occurrences;

//     //  Create list of k-mers with occurrence numbers
//     for(size_t i = 0; i < seq.size() - k + 1; ++i) {
//         auto occ = occurrences[seq.substr(i, k)]++;
//         mers.emplace_back(i, occ, (uint64_t)0);
//     }

//     xxhash hash;
//     std::vector<uint64_t> omh_results;
//     for(unsigned i = 0; i < m; ++i) {
//         const auto seed = prg();
//         for(auto& meri : mers) {
//         hash.reset(seed);
//         hash.update(&seq.data()[meri.pos], k);
//         if(weight) hash.update(&meri.occ, sizeof(meri.occ));
//         meri.hash = hash.digest();
//         }

//         // std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), & { return x.hash < y.hash; });
//         // std::partial_sort(mers.begin(), mers.begin() + l, mers.end(),  { return x.hash < y.hash; });

//         std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
//         std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });


//         omh_results.push_back(mers[0].hash); // Save the smallest hash as the OMH result for this round
//     }

//     return omh_results; // Return the OMH results
// }


#endif /* __OMH_H__ */
