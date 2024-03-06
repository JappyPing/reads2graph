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
// #include <seqan3/utility/views/pairwise_combine.hpp>
#include <boost/functional/hash.hpp>
// #include <xxhash.h>
// #include <seqan3/alphabet/nucleotide/dna5.hpp>
// #include <seqan3/alphabet/range/hash.hpp>
// #include <seqan3/utility/container/dynamic_bitset.hpp>
// #include <seqan3/search/views/kmer_hash.hpp>
// #include <seqan3/core/debug_stream.hpp>
// #include <seqan3/alphabet/hash.hpp>
// #include <seqan3/alphabet/hash.hpp>
// #include <seqan3/alphabet/range/hash.hpp>
// OMH::OMH(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : read2count(read2count), args(args) {}
OMH::OMH(cmd_arguments args) : args(args) {}

unsigned OMH::omh_k(unsigned L, double p, uint8_t d) {
    unsigned k = ceil((p*(1+L))/(d+p));
    if (k < 4){
        k = 4;
        Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Better k {} has been changed to 4.", k));
    } else if (k > (args.read_length / 4)) {
        k = args.read_length/4;
        Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Better k {} has been changed to {}.", k, args.read_length/4));                
    }
    return k;
}

std::vector<std::pair<std::uint64_t, unsigned>> OMH::get_seeds_k(){
    unsigned k;
    if (args.minimizer_omh) {
        k = args.omh_k + args.omh_k_step_size;
    } else {
        k = omh_k(args.read_length, args.bad_kmer_ratio/3, 1);
    } 
    // Use a random_device to seed the random number generator
    std::random_device rd;
    // Use the Mersenne Twister engine for random number generation
    std::mt19937_64 generator(rd());
    // Specify the range of values for your seeds
    std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());

    // Generate seeds and store them in a vector
    std::vector<std::pair<std::uint64_t, unsigned>> seeds_k;
    for (unsigned int i = 0; i < args.omh_times; ++i) {
        auto cur_seed = distribution(generator);
        Utils::getInstance().logger(LOG_LEVEL_INFO, std::format("k: {} and seed: {};", k, cur_seed));
        std::pair<std::uint64_t, unsigned> cur_pair = std::make_pair(cur_seed, k);
        seeds_k.push_back(cur_pair);
        // k = k + args.omh_k_step_size;
        k++;
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, "The above k and seed pairs are used for OMH bucketing.");
    return seeds_k;
}

std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH::omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k){
    // int available_cores = omp_get_max_threads();
    // auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    // omp_set_num_threads(num_cores_to_use);
    #pragma omp parallel for
    for (auto const & read : unique_reads){
        
        for(auto &pair : seeds_k){
            std::uint64_t seed = pair.first;
            unsigned k = pair.second;
            auto omh_value = omh_pos(read, k, seed); 
            #pragma omp critical
            {
                omh2reads[omh_value].push_back(read);
            }        
        } 
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
    return omh2reads;            
}

// uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, std::uint64_t seed) {
uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed) {
    if(read.size() < k) return {};
    std::vector<std::uint64_t> hash_vec;
    std::unordered_map<std::string, unsigned> occurrences;
    std::uint64_t cur_seed = seed;

    auto seql = read | seqan3::views::to_char;
    string read_str(seql.begin(), seql.end());
    for(size_t i = 0; i < read_str.size() - k + 1; ++i) {
        string kmer = read_str.substr(i, k);
        occurrences[kmer]++;
        boost::hash_combine(cur_seed, kmer);
        boost::hash_combine(cur_seed, occurrences[kmer]);
        hash_vec.emplace_back(cur_seed);
        cur_seed = seed;
    }
    auto min_hash = std::min_element(hash_vec.begin(), hash_vec.end());    
    return *min_hash;
}


// std::vector<std::pair<std::uint64_t, unsigned>> OMH::get_seeds_k(){
//     // Use a random_device to seed the random number generator
//     std::random_device rd;
//     // Use the Mersenne Twister engine for random number generation
//     std::mt19937_64 generator(rd());
//     // Specify the range of values for your seeds
//     std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());

//     // Generate seeds and store them in a vector
//     std::vector<std::pair<std::uint64_t, unsigned>> seeds_k;
//     auto dis = args.max_edit_dis;
//     unsigned omh_times;
//     if (args.max_edit_dis > 2){
//         omh_times = 1;
//     } else {
//         omh_times = 2;
//     }
//     for (unsigned int i = 0; i < args.max_edit_dis; ++i) {
//         for (unsigned int j = 0; j < omh_times; ++j){
//             auto cur_seed = distribution(generator); 
//             auto cur_k = omh_k(args.read_length, args.bad_kmer_ratio/3, dis);
//             Utils::getInstance().logger(LOG_LEVEL_INFO, std::format("k={} and seed={} are used for OMH bucketing.", cur_k, cur_seed));
//             std::pair<std::uint64_t, unsigned> cur_pair = std::make_pair(cur_seed, cur_k);
//             seeds_k.push_back(cur_pair);
//         }
//         dis--;
//     }

//     // unsigned cur_k = betterK;
//     // unsigned even_k = betterK;
//     // unsigned odd_k = betterK;
//     // auto dis = args.max_edit_dis;
//     // for (unsigned int i = 0; i < args.omh_times; ++i) {
//     //     auto cur_seed = distribution(generator);
//     //     unsigned cur_k;
//     //     if (dis > 0){
//     //         cur_k = omh_k(args.read_length, args.bad_kmer_ratio, dis);
//     //         dis--;
//     //     }
//     //     if (dis == 0 && i < args.omh_times){
//     //         cur_k = omh_k(args.read_length, args.bad_kmer_ratio, 1);
//     //         dis--;
//     //     }

//     //     Utils::getInstance().logger(LOG_LEVEL_INFO, std::format("k={} and seed={} are used for OMH bucketing.", cur_k, cur_seed));
//     //     std::pair<std::uint64_t, unsigned> cur_pair = std::make_pair(cur_seed, cur_k);
//     //     seeds_k.push_back(cur_pair);
//     //     // if (cur_k >= 4 && cur_k < (args.read_length - 2)) {
//     //     //     if (i % 2 == 0){
//     //     //         cur_k = even_k + args.omh_k_step_size;
//     //     //         even_k = cur_k;
//     //     //     } else {
//     //     //         cur_k = odd_k - args.omh_k_step_size;
//     //     //         odd_k = cur_k;
//     //     //     }
//     //     // } else {
//     //     //     cur_k = betterK;
//     //     // }
//     // }
//     Utils::getInstance().logger(LOG_LEVEL_DEBUG,  "Seeds and kmer-size done!");
//     return seeds_k;
// }


// std::vector<std::pair<std::uint64_t, unsigned>> OMH::get_seeds_k(){
//     auto betterK = omh_k(args.read_length, args.bad_kmer_ratio, args.max_edit_dis);

//     Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Better k for bucketing by OMH: {}.", betterK));  
//     // std::default_random_engine prg;
//     // vector<unsigned int> seeds;
//     // for(unsigned i = 0; i < args.omh_times; ++i) {
//     //     seeds.push_back(prg());
//     // }
//     // Use a random_device to seed the random number generator
//     std::random_device rd;
//     // Use the Mersenne Twister engine for random number generation
//     std::mt19937_64 generator(rd());
//     // Specify the range of values for your seeds
//     std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());

//     // Generate seeds and store them in a vector
//     std::vector<std::pair<std::uint64_t, unsigned>> seeds_k;
//     unsigned cur_k = betterK;
//     unsigned even_k = betterK;
//     unsigned odd_k = betterK;
//     for (unsigned int i = 0; i < args.omh_times; ++i) {
//         auto cur_seed = distribution(generator);
//         Utils::getInstance().logger(LOG_LEVEL_INFO, std::format("k={} and seed={} are used for OMH bucketing.", cur_k, cur_seed));
//         std::pair<std::uint64_t, unsigned> cur_pair = std::make_pair(cur_seed, cur_k);
//         seeds_k.push_back(cur_pair);
//         if (cur_k >= 4 && cur_k < (args.read_length - 2)) {
//             if (i % 2 == 0){
//                 cur_k = even_k + args.omh_k_step_size;
//                 even_k = cur_k;
//             } else {
//                 cur_k = odd_k - args.omh_k_step_size;
//                 odd_k = cur_k;
//             }
//         } else {
//             cur_k = betterK;
//         }
//     }

//     // for(auto &pair : seeds_k){
//     //     std::uint64_t seed = pair.first;
//     //     unsigned k = pair.second;
//     //     Utils::getInstance().logger(LOG_LEVEL_INFO, std::format("k={} and seed={} are used for OMH bucketing.", k, seed));
//     // }
//     Utils::getInstance().logger(LOG_LEVEL_DEBUG,  "Seeds and kmer-size done!");
//     return seeds_k;
// }

// std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH::omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k){
//     int available_cores = omp_get_max_threads();
//     auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
//     omp_set_num_threads(num_cores_to_use);

//     // auto uniq_num = read2count.size();
//     // auto uniq_num = unique_reads.size();
    
//     // #pragma omp parallel for
//     // for (size_t i = 0; i < uniq_num; ++i) {
//     //     auto it = std::next(read2count.begin(), i);
//     //     const auto& [read, count] = *it;
//     #pragma omp parallel for
//     for (auto const & read : unique_reads){
//         // std::vector<uint64_t> cur_omh_vals;
        
//         for(auto &pair : seeds_k){
//             std::uint64_t seed = pair.first;
//             unsigned k = pair.second;
//             // auto omh_value = omh_pos(read, 28, args.omh_kmer_n, seed);   
//             auto omh_value = omh_pos(read, k, args.omh_kmer_n, seed);  
//             // auto omh_value = omh_pos(read, args.omh_k, seed);  
//             // cur_omh_vals.push_back(omh_value); 
//             auto omh_min_max = omh_pos(read, k, args.omh_kmer_n, seed); 
//             #pragma omp critical
//             {
//                 omh2reads[omh_value].push_back(read);
//                 // omh2reads[omh_min_max.first].push_back(read);
//                 // omh2reads[omh_min_max.second].push_back(read);
//             }        
//         } 
//     }
//         // auto omh_combins = seqan3::views::pairwise_combine(cur_omh_vals);
//         // for (size_t i = 0; i < omh_combins.size(); ++i)
//         // {
//         //     auto const &omh_comb = omh_combins[i];
//         //     auto const &val1 = std::get<0>(omh_comb);
//         //     auto const &val2 = std::get<1>(omh_comb);
//         //     std::size_t new_val = 0;
//         //     boost::hash_combine(new_val, val1);
//         //     boost::hash_combine(new_val, val2);
//         //     // auto new_val = val1 ^ val2;
//         //     #pragma omp critical
//         //     {
//         //         omh2reads[new_val].push_back(read);
//         //     
//         // }
//     // }
//     Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
//     return omh2reads;            
// }

// // uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned int seed) {
// // uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, std::uint64_t seed) {
// uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, std::uint64_t seed) {

//     if(read.size() < k) return {};

//     const bool weight = l > 0;
//     // const bool weight = 1 > 0;
//     if(l == 0) l = 1;

//     std::vector<mer_info> mers;
//     std::vector<std::uint64_t> hash_vec;
//     // std::unordered_map<std::uint64_t, unsigned> occurrences;
//     std::unordered_map<std::string, unsigned> occurrences;
//     std::uint64_t cur_seed = seed;
//     // in order to hash a kmer using seqan3::views::kmer_hash, here I use k+1 rather than k to calculate the kmers for a reads. using k the shape size may be the same as the kmer size which will be terminated by seqan3

//     auto seql = read | seqan3::views::to_char;
//     string read_str(seql.begin(), seql.end());
//     for(size_t i = 0; i < read_str.size() - k + 1; ++i) {
//         string kmer = read_str.substr(i, k);
//         // auto occ = occurrences[kmer]++;
//         occurrences[kmer]++;
//         boost::hash_combine(cur_seed, kmer);
//         if (weight)
//             boost::hash_combine(cur_seed, occurrences[kmer]);
//         // mers.emplace_back(i, occ, cur_seed);
//         hash_vec.emplace_back(cur_seed);
//         cur_seed = seed;
//     }

//     // for (size_t i = 0; i < read.size() - k + 1; ++i)
//     // {
//     //     auto kmer = std::vector<seqan3::dna5>(read.begin() + i, read.begin() + i + k);
//     //     auto hash_val_range = kmer | seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(k)});
//     //     std::uint64_t hash_val = *hash_val_range.begin();

//     //     auto occ = occurrences[hash_val]++;
//     //     // seqan3::debug_stream << kmer.size() << " " << k << " " << hash_val_range.size() << '\n';        
//     //     boost::hash_combine(cur_seed, hash_val);
//     //     if (weight)
//     //         boost::hash_combine(cur_seed, occurrences[hash_val]);
//     //     // seqan3::debug_stream << i << " " << kmer << " " << hash_val << " " << cur_seed << " " << occ << " " << cur_seed << '\n';
//     //     mers.emplace_back(i, occ, cur_seed);
//     //     // minHash = std::min(minHash, cur_seed);
//     //     cur_seed = seed;
//     // }
//     // return minHash;
//     // std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
//     // std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });

//     // return mers[0].hash; // Return the OMH results
//     auto min_hash = std::min_element(hash_vec.begin(), hash_vec.end());    
//     return *min_hash;
//     // auto minMaxPair = std::minmax_element(hash_vec.begin(), hash_vec.end());
//     // std::pair<uint64_t, uint64_t> min_max = std::make_pair(*minMaxPair.first, *minMaxPair.second);
//     // return min_max;
// }

// uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned l, unsigned int seed) {
// // uint64_t OMH::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, unsigned int seed) {
//     auto seql = read | seqan3::views::to_char;
//     string seq(seql.begin(), seql.end());
//     if(seq.size() < k) return {};

//     const bool weight = l > 0;
//     // const bool weight = 1 > 0;
//     if(l == 0) l = 1;

//     std::vector<mer_info> mers;
//     std::unordered_map<std::string, unsigned> occurrences;

//     //  Create list of k-mers with occurrence numbers
//     for(size_t i = 0; i < seq.size() - k + 1; ++i) {
//         auto occ = occurrences[seq.substr(i, k)]++;
//         mers.emplace_back(i, occ, (uint64_t)0);
//     }

//     xxhash hash;
//     // std::uint64_t minHash = std::numeric_limits<uint64_t>::max();  // Initialize with the maximum value
//     for(auto& meri : mers) {
//         hash.reset(seed);
//         hash.update(&seq.data()[meri.pos], k);
//         if(weight) hash.update(&meri.occ, sizeof(meri.occ));
//         meri.hash= hash.digest();
//         // minHash = std::min(minHash, meri.hash); // using minhash directly will get poor performance
//     }
//     // return minHash;
//     std::partial_sort(mers.begin(), mers.begin() + l, mers.end(), [&](const mer_info& x, const mer_info& y) { return x.hash < y.hash; });
//     std::sort(mers.begin(), mers.begin() + l, [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });

//     return mers[0].hash; // Return the OMH results
// }
// std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH::omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k){
//     int available_cores = omp_get_max_threads();
//     auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
//     omp_set_num_threads(num_cores_to_use);

//     // auto uniq_num = read2count.size();
//     // auto uniq_num = unique_reads.size();
    
//     // #pragma omp parallel for
//     // for (size_t i = 0; i < uniq_num; ++i) {
//     //     auto it = std::next(read2count.begin(), i);
//     //     const auto& [read, count] = *it;
//     #pragma omp parallel for
//     for (auto const & read : unique_reads){
//         // std::vector<uint64_t> cur_omh_vals;
        
//         for(auto &pair : seeds_k){
//             std::uint64_t seed = pair.first;
//             unsigned k = pair.second;
//             // auto omh_value = omh_pos(read, 28, args.omh_kmer_n, seed);   
//             // auto omh_value = omh_pos(read, k, args.omh_kmer_n, seed);  
//             // auto omh_value = omh_pos(read, args.omh_k, seed);  
//             // cur_omh_vals.push_back(omh_value); 
//             auto omh_min_max = omh_pos(read, k, seed); 
//             #pragma omp critical
//             {
//                 omh2reads[omh_value].push_back(read);
//                 // omh2reads[omh_min_max.first].push_back(read);
//                 // omh2reads[omh_min_max.second].push_back(read);
//             }        
//         } 
//     }
//         // auto omh_combins = seqan3::views::pairwise_combine(cur_omh_vals);
//         // for (size_t i = 0; i < omh_combins.size(); ++i)
//         // {
//         //     auto const &omh_comb = omh_combins[i];
//         //     auto const &val1 = std::get<0>(omh_comb);
//         //     auto const &val2 = std::get<1>(omh_comb);
//         //     std::size_t new_val = 0;
//         //     boost::hash_combine(new_val, val1);
//         //     boost::hash_combine(new_val, val2);
//         //     // auto new_val = val1 ^ val2;
//         //     #pragma omp critical
//         //     {
//         //         omh2reads[new_val].push_back(read);
//         //     
//         // }
//     // }
//     Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
//     return omh2reads;            
// }