// OMH.hpp

/*The implementation of Order Min Hash (OMH) is modified based on the original implementation of OMH(https://github.com/Kingsford-Group/omhismb2019). If you want to use the relevant source codes in your project, please remember to cite the original work listed below.
Guillaume Marçais, Dan DeBlasio, Prashant Pandey, Carl Kingsford, Locality-sensitive hashing for the edit distance, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i127–i135, https://doi.org/10.1093/bioinformatics/btz354
*/

#ifndef __GOMH_HPP__
#define __GOMH_HPP__

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

class gOMH
{
public:
    // OMH(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    gOMH(cmd_arguments args);
    // std::vector<std::pair<std::uint64_t, unsigned>> get_seeds_k();
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k);
    uint64_t gomh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed);
    
    unsigned gomh_k(unsigned L, double p, uint8_t d);
    std::string getGappedSubstring(const std::string& str, size_t startPos, size_t length);
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::uint64_t seed, unsigned k);
    std::vector<std::uint64_t> gomh_pos2(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed);
    ////////////////////////////////
    // std::vector<uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed, int m);
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omhs2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned k, std::uint64_t seed, int m);
    ////////////////////////////////
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::uint64_t> seeds, unsigned k);
    
    // uint64_t omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned int seed);
private:
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2reads;
    cmd_arguments args;

};

/////////////////////////////////////////////////////////////////////////////////////////
// This is the source code for the function of omh computation from the authors, please cite the following the literature if you want to use the source code
// Reference: Guillaume Marçais, Dan DeBlasio, Prashant Pandey, Carl Kingsford, Locality-sensitive hashing for the edit distance, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i127–i135, https://doi.org/10.1093/bioinformatics/btz354
// Compute the position in sequence of the k-mers picked by omh, and
// passed them 1 by 1 to block. block takes 3 arguments: i \in [m], j
// \in [l] and the position.
// template<typename EngineT, typename BT>
// void omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, EngineT& prg, BT block) {
template<typename EngineT>
std::vector<uint64_t> ori_omh_pos(const std::string& seq, unsigned k, unsigned l, unsigned m, EngineT& prg) {
    
  if(seq.size() < k) return;
  const bool weight = l > 0;
  if(l == 0) l = 1;

  std::vector<mer_info> mers;
  std::unordered_map<std::string, unsigned> occurrences;
  size_t pos[l];

  //  Create list of k-mers with occurrence numbers
  for(size_t i = 0; i < seq.size() - k + 1; ++i) {
    auto occ = occurrences[seq.substr(i, k)]++;
    mers.emplace_back(i, occ, (uint64_t)0);
  }

  xxhash hash;
  for(unsigned i = 0; i < m; ++i) {
    const auto seed = prg();
    for(auto& meri : mers) {
      hash.reset(seed);
      hash.update(&seq.data()[meri.pos], k);
      if(weight) hash.update(&meri.occ, sizeof(meri.occ));
      meri.hash = hash.digest();
    }

    std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
    std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });
//     for(unsigned j = 0; j < l; ++j)
//       block(i, j, mers[j].pos);
    }
    uint64_t smallest_hash = mers[0].hash;
    return smallest_hash;
}


#endif /* __OMH_H__ */
