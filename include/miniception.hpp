/*
The head is used to implementate the miniception c++ version based on the source codes at https://github.com/myprogrammerpersonality/Minimizers/tree/main. If you want to use the following source codes in your project, please remember to cite the original work (see below) and the implementation of miniception c++ version at https://github.com/myprogrammerpersonality/Minimizers/tree/main.
Hongyu Zheng, Carl Kingsford, Guillaume Marçais, Improved design and analysis of practical minimizers, Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i119–i127, https://doi.org/10.1093/bioinformatics/btaa472
*/
#ifndef __MINICEPTION_H__
#define __MINICEPTION_H__

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <xxhash.h>
#include <seeded_prg.hpp>
#include <set>
#include <deque>
#include <utility>
#include <functional>
#include <seqan3/std/ranges>

// MonotonicQueue implementation
template <typename X, typename T>
class MonotonicQueue {
private:
    std::deque<std::pair<X, T>> dq;

public:
    void Insert(X item, T time) {
        while (!dq.empty() && dq.back().first > item) {
            dq.pop_back();
        }
        dq.emplace_back(item, time);
    }

    std::pair<X, T> Fetch(T time) {
        while (!dq.empty() && dq.front().second < time) {
            dq.pop_front();
        }
        return dq.empty() ? std::pair<X, T>() : dq.front();
    }
};

class Miniception {
public:
    // Function to group reads by minimizers using a sliding window
    template<typename EngineT>
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> miniception2read_main(
        std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned kmer_size, unsigned window_size, EngineT& prg, unsigned num_process){
        std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> miniception2reads;

        const auto seed = prg();
        xxhash hash; // Initialize xxHash

        #pragma omp parallel for num_threads(num_process) schedule(static)
        for (size_t idx = 0; idx < unique_reads.size(); ++idx) {
            const auto& read = unique_reads[idx];
            auto seql = read | seqan3::views::to_char;
            std::string read_str(seql.begin(), seql.end());

            MonotonicQueue<std::uint64_t, int> mq;
            std::vector<std::uint64_t> read_minimizers;

            // Compute minimizers for the current read
            for (size_t i = 0; i <= read_str.size() - kmer_size; ++i) {
                std::string kmer = read_str.substr(i, kmer_size);
                
                // Compute xxHash for the k-mer
                hash.reset(seed); 
                hash.update(kmer.data(), kmer.size());
                std::uint64_t hash_val = hash.digest();

                mq.Insert(hash_val, i);

                if (i >= window_size - 1) {
                    auto fetch_pair = mq.Fetch(i - window_size + 1);
                    if (read_minimizers.empty() || read_minimizers.back() != fetch_pair.first) {
                        read_minimizers.push_back(fetch_pair.first);
                    }
                }
            }

            // Group reads by their minimizer hashes
            #pragma omp critical
            {
                for (auto min_hash : read_minimizers) {
                    miniception2reads[min_hash].push_back(read);
                }
            }
        }

        return miniception2reads;
    }

};

#endif // __MINICEPTION_H__
