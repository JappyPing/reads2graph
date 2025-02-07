// OMH.hpp

/*The implementation of Order Min Hash (OMH) is modified based on the original implementation of OMH(https://github.com/Kingsford-Group/omhismb2019). If you want to use the relevant source codes in your project, please remember to cite the original work listed below.
Guillaume Marçais, Dan DeBlasio, Prashant Pandey, Carl Kingsford, Locality-sensitive hashing for the edit distance, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i127–i135, https://doi.org/10.1093/bioinformatics/btz354
*/

#ifndef __GOMH_HPP__
#define __GOMH_HPP__

#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include "omh.hpp"

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

// struct mer_info {
//   size_t pos;
//   std::uint64_t hash;
//   unsigned occ;
//   mer_info(size_t p, unsigned o, std::uint64_t h)
//     : pos(p)
//     , hash(h)
//     , occ(o)
//   { }
// };

class gOMH
{
public:
    // OMH(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    gOMH(cmd_arguments args);
    // std::vector<std::pair<std::uint64_t, unsigned>> get_seeds_k();
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k);
    std::uint64_t gomh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed);
    
    unsigned gomh_k_size(unsigned L, double p, uint8_t d);
    std::string getGappedSubstring(const std::string& str, size_t startPos, size_t length);
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::uint64_t seed, unsigned k);
    std::vector<std::uint64_t> gomh_pos2(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed);
    ////////////////////////////////
    // std::vector<std::uint64_t> omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed, int m);
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omhs2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned k, std::uint64_t seed, int m);
    ////////////////////////////////
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::uint64_t> seeds, unsigned k);
    
    // std::uint64_t omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned int seed);
private:
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2reads;
    cmd_arguments args;

};
#endif /* __OMH_H__ */
