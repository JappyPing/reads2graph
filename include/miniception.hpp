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
#include <seqan3/alphabet/all.hpp>
#include <xxhash.h>
#include <seeded_prg.hpp>
#include <set>
#include <deque>
#include <utility>
#include <functional>

// MonotonicQueue implementation at https://github.com/myprogrammerpersonality/Minimizers/tree/main
template <typename X, typename T>
class MonotonicQueue {
    private:
        std::deque<std::pair<X, T>> dq;

    public:
        // Insert an element into the queue
        void Insert(X item, T time) {
            while (!dq.empty() && dq.back().first > item) {
                dq.pop_back();
            }
            dq.emplace_back(item, time);
        }

        // Fetch the current minimum and remove old elements
        std::pair<X, T> Fetch(T time) {
            while (!dq.empty() && dq.front().second < time) {
                dq.pop_front();
            }
            return dq.empty() ? std::pair<X, T>() : dq.front();
        }
};

class Miniception {
public:
    template<typename EngineT>
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> miniception2read_main(
        std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned kmer_size, unsigned window_size, EngineT& prg, unsigned num_process){

        std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> miniception2reads;
        std::cout << unique_reads.size() << std::endl;
        const auto seed = prg();
         
        // std::hash<std::string> str_hash;
        #pragma omp parallel for num_threads(num_process) schedule(static)
        for (const auto& read : unique_reads) {      
            auto seql = read | seqan3::views::to_char;
            // string read_str(seql.begin(), seql.end());
            std::vector<char> read_vec(seql.begin(), seql.end());

            MonotonicQueue<std::uint64_t, int> mq;
            std::set<std::uint64_t> minimizer_val;
            xxhash hash;
            // Compute minimizers for the current read
            // for (size_t i = 0; i <= read_str.size() - kmer_size; ++i){
            for (size_t i = 0; i <= read_vec.size() - kmer_size; ++i) {
                // std::string kmer = read_str.substr(i, kmer_size);
                std::vector<char> kmer(read_vec.begin() + i, read_vec.begin() + i + kmer_size);
                
                // Compute xxHash for the k-mer
                hash.reset(seed); 
                hash.update(kmer.data(), kmer.size());
                std::uint64_t hash_val = hash.digest();
                // size_t hash_val = static_cast<uint64_t>(str_hash(kmer));
                mq.Insert(hash_val, i);

                if (i >= window_size - 1) {
                    auto fetch_pair = mq.Fetch(i - window_size + 1);
                    minimizer_val.insert(fetch_pair.first);    
                }
            }
            std::vector<std::uint64_t> minimizers(minimizer_val.begin(), minimizer_val.end());
            // Group reads by their minimizer hashes
            for (auto min_hash : minimizers) {
                #pragma omp critical
                {
                    miniception2reads[min_hash].push_back(read);
                }  
            }
        }
        std::cout << "test" << miniception2reads.size() << std::endl;
        return miniception2reads;
    }

};

#endif // __MINICEPTION_H__
