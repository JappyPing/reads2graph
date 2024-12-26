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
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> miniception2read_main(
        std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned kmer_size, unsigned window_size, unsigned num_process);

    std::vector<std::uint64_t> miniception_main(std::vector<seqan3::dna5> read, unsigned kmer_size, unsigned window_size, uint64_t seed);
};

#endif // __MINICEPTION_H__
