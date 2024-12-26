#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include "miniception.hpp"

Miniception::Miniception(cmd_arguments args) : args(args) {}

std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> Miniception::miniception2read_main(
    std::vector<std::vector<seqan3::dna5>> unique_reads, unsigned kmer_size, unsigned window_size, unsigned num_process){

    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> miniception2reads;
    std::cout << unique_reads.size() << std::endl;

    // std::hash<std::string> str_hash;
    #pragma omp parallel for num_threads(num_process) schedule(static)
    for (const auto& read : unique_reads) {      
        auto minimizers = miniception_main(read, kmer_size, window_size, args.seed);
        // Group reads by their minimizer hashes
        for (auto min_hash : minimizers) {
            #pragma omp critical
            {
                miniception2reads[min_hash].push_back(read);
            }  
        }
    }
    //std::cout << "test" << miniception2reads.size() << std::endl;
    return miniception2reads;
}

std::vector<std::uint64_t> Miniception::miniception_main(std::vector<seqan3::dna5> read, unsigned kmer_size, unsigned window_size, uint64_t seed){
    auto seql = read | seqan3::views::to_char;
    std::vector<char> read_vec(seql.begin(), seql.end());

    MonotonicQueue<std::uint64_t, int> mq;
    std::set<std::uint64_t> minimizer_val;
    xxhash hash;
    for (size_t i = 0; i <= read_vec.size() - kmer_size; ++i) {
        std::vector<char> kmer(read_vec.begin() + i, read_vec.begin() + i + kmer_size);
        
        // Compute xxHash for the k-mer
        hash.reset(seed); 
        hash.update(kmer.data(), kmer.size());
        std::uint64_t hash_val = hash.digest();
        mq.Insert(hash_val, i);

        if (i >= window_size - 1) {
            auto fetch_pair = mq.Fetch(i - window_size + 1);
            minimizer_val.insert(fetch_pair.first);    
        }
    }
    std::vector<std::uint64_t> minimizers(minimizer_val.begin(), minimizer_val.end());
    return minimizers;
}