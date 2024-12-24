/*
The head is used to implementate the Order Min Hash (OMH)(https://github.com/Kingsford-Group/omhismb2019). If you want to use the following source codes in your project, please remember to cite the original work listed below.
Guillaume Marçais, Dan DeBlasio, Prashant Pandey, Carl Kingsford, Locality-sensitive hashing for the edit distance, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i127–i135, https://doi.org/10.1093/bioinformatics/btz354
*/

#ifndef __OMH_H__
#define __OMH_H__

#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <limits>
#include <seqan3/alphabet/all.hpp>
#include <xxhash.hpp>
#include <seeded_prg.hpp>
#include <omp.h>

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

// Compute the position in sequence of the k-mers picked by omh, and
// passed them 1 by 1 to block. block takes 3 arguments: i \in [m], j
// \in [l] and the position.
template<typename EngineT, typename BT>
void omh_pos(const std::string& seq, unsigned k, unsigned l, EngineT& prg, BT block) {
    if (seq.size() < k) return;
    // const bool weight = l > 0;
    if (l == 0) l = 1;

    std::vector<mer_info> mers;
    std::unordered_map<std::string, unsigned> occurrences;

    // Create list of k-mers with occurrence numbers
    for (size_t i = 0; i < seq.size() - k + 1; ++i) {
        auto occ = occurrences[seq.substr(i, k)]++;
        mers.emplace_back(i, occ, (uint64_t)0);
    }

    // Sort k-mers by position
    std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), 
        [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });

    // Call the block function on the sorted k-mers
    for (unsigned j = 0; j < l; ++j) {
        block(mers[j].pos);
    }
}

struct sketch {
    std::string name;       // Optional name for the sketch
    unsigned k, l, m;       // Parameters for sketching
    std::vector<char> data; // Data for the sketch (k-mers)
};

template<typename EngineT = std::mt19937_64>
class omh_sketcher {
protected:
    unsigned m_k, m_l;         // Parameters for sketching
    uint64_t m_seed;           // Seed for PRG
    EngineT m_prg;             // Random number generator (PRG)

    // Compute sketch positions for the given sequence and store them in 'ptr'
    inline void compute_sketch(uint8_t* ptr, const char* seq) {
        omh_pos(seq, m_k, m_l, m_prg, [&ptr, &seq, this](size_t pos) {
            std::memcpy(ptr, seq + pos, m_k);  // Copy k-mer into sketch
            ptr += m_k;
        });
    }

public:
    // Constructor: Initialize k, l, values and the PRG with a fixed seed
    omh_sketcher(unsigned k, unsigned l, uint64_t seed)
        : m_k(k), m_l(l), m_seed(seed), m_prg(seed) // Pass seed to PRG
    { }

    // Compute sketch from sequence
    void compute(const std::string& seq, sketch& sk) {
        sk.k = m_k;
        sk.l = m_l;
        sk.data.resize(std::max(sk.l, (unsigned)1) * sk.k);

        // Use the same seed for each read to keep sketch consistent
        m_prg.seed(m_seed);  
        compute_sketch(reinterpret_cast<uint8_t*>(sk.data.data()), seq.data());
    }

    // Compute and return sketch directly from sequence
    sketch compute(const std::string& seq) {
        sketch sk;
        compute(seq, sk);
        return sk;
    }
};

class OMH {
public:
    // Function to group reads by OMH sketch using m hash functions
    template<typename EngineT>
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>>
    omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, omh_sketcher<EngineT>& sketcher, EngineT& prg, unsigned m, int num_process) {
        std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> omh2reads;
        xxhash hash;

        // Iterate over m hash functions (or seeds)
        for (unsigned i = 0; i < m; ++i) {
            // Generate a new seed for this hash iteration
            const auto seed = prg();

            #pragma omp parallel for num_threads(num_process) schedule(static)
            for (const auto& read : unique_reads) {

                auto seql = read | seqan3::views::to_char;
                string read_str(seql.begin(), seql.end());

                // Compute the sketch for the current read
                sketch sk = sketcher.compute(read_str);

                // Hash the sketch to generate a uint64_t key using the current seed
                hash.reset(seed);
                hash.update(sk.data.data(), sk.data.size());
                std::uint64_t sketch_hash = hash.digest();

                // Group reads by their sketch hash
                #pragma omp critical
                {
                    omh2reads[sketch_hash].push_back(read);
                }
            }
        }

        return omh2reads;
    }
};


#endif /* __OMH_H__ */