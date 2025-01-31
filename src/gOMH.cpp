#include "gOMH.hpp"
#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
// #include <format>
#include <boost/format.hpp>
#include <seqan3/alphabet/all.hpp>
// #include <seqan3/utility/views/pairwise_combine.hpp>
#include <boost/functional/hash.hpp>

gOMH::gOMH(cmd_arguments args) : args(args) {}

unsigned gOMH::gomh_k(unsigned L, double p, uint8_t d) {
    unsigned k;
    if (args.default_params) {
        if (args.read_length >= 6 && args.read_length < 16){ 
            k = 3;
        } else {
            k = ceil(((1-p)*(2+L))/(d+2-2*p));
            if (k < 4){
                k = 4;
            } else if (k > 27) {
                k = 27;               
            }      
            auto gomh_kmer_n = L - 2 * k + 1;
            while ((gomh_kmer_n <= args.gomh_times + 1) && k > 4){
                k--;
            }
        }

    } else {
        k = args.gomh_k;
    }        
    return k;
}
/*
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
        if (args.read_length > 50) {
           k++; 
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, "The above k and seed pairs are used for OMH bucketing.");
    return seeds_k;
}
*/
std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gOMH::gomh2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::uint64_t seed, unsigned k){
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto const & read : unique_reads){
        auto gomh_value = gomh_pos(read, k, seed); 
        #pragma omp critical
        {
            gomh2reads[gomh_value].push_back(read);
        }        
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for gOMH done.");  
    return gomh2reads;            
}

// std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH::omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::uint64_t> seeds, unsigned k){

//     #pragma omp parallel for num_threads(args.num_process) schedule(static)
//     for (auto const & read : unique_reads){
//         #pragma omp parallel for num_threads(args.num_process) schedule(static)
//         for(auto &seed : seeds){
//             auto omh_value = omh_pos(read, k, seed); 
//             #pragma omp critical
//             {
//                 omh2reads[omh_value].push_back(read);
//             }        
//         } 
//     }
//     Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
//     return omh2reads;            
// }


std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gOMH::gomh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k){
    std::size_t reads_n = unique_reads.size();
    if (args.gomh_flag){
        // When the permutation_times larger than the number of k-mer candidates and the kmer size are the same one, bucketing the reads using each kmer candidate
        auto first_pair = seeds_k[0];
        std::uint64_t seed = first_pair.first;
        unsigned k = first_pair.second;
        while (1) {
            std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2reads;
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (auto const & read : unique_reads){  
                auto gomh_values = gomh_pos2(read, k, seed);
                for (auto const & gomh_val : gomh_values){
                    #pragma omp critical
                    {
                        gomh2reads[gomh_val].push_back(read);
                    }                    
                }
            }
            std::size_t bin_n = gomh2reads.size();
            if (bin_n < (reads_n * args.bin2reads_ratio)){
                return gomh2reads;  
            } else {
                k -= 2;
                if (k == 4){
                    return gomh2reads;
                }
            }                
        } 
    } else {
        std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> gomh2reads;
        for(auto &pair : seeds_k){
            std::uint64_t seed = pair.first;
            unsigned k = pair.second;
            while (1) {
                std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> temp_gomh2reads;
                #pragma omp parallel for num_threads(args.num_process) schedule(static)
                for (auto const & read : unique_reads){

                    auto gomh_value = gomh_pos(read, k, seed); 
                    #pragma omp critical
                    {
                        temp_gomh2reads[gomh_value].push_back(read);
                    }        
                } 
                std::size_t bin_n = temp_gomh2reads.size();
                if (bin_n < (reads_n * args.bin2reads_ratio)){
                    for (auto& [key, value] : temp_gomh2reads) {
                        gomh2reads[key].insert(gomh2reads[key].end(), 
                                            std::make_move_iterator(value.begin()), 
                                            std::make_move_iterator(value.end()));
                    }
                    temp_gomh2reads.clear();              
                    break;  
                } else {
                    k -= 2;
                    if (k == 4){
                        for (auto& [key, value] : temp_gomh2reads) {
                            gomh2reads[key].insert(gomh2reads[key].end(), 
                                                std::make_move_iterator(value.begin()), 
                                                std::make_move_iterator(value.end()));
                        }
                        temp_gomh2reads.clear();                             
                        break;
                    }
                }                       
            }
        } 
        Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for gOMH done.");
        return gomh2reads;       
    }
    // std::size_t bin_n = gomh2reads.size();
    // if (bin_n < (reads_n * args.bin2reads_ratio)){
    //     Utils::getInstance().logger(LOG_LEVEL_DEBUG, boost::str(boost::format("Size of gomh2reads: %1%!") % bin_n));
    //     return minimiser2reads;  
    // } else {
    //     k -= 2;
    //     if (k == 4){
    //         return minimiser2reads;
    //     }
    // }

    // Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for gOMH done.");  
    // return gomh2reads;
}            

std::string gOMH::getGappedSubstring(const std::string& str, size_t startPos, size_t length) {
    std::string gappedSubstring;

    for (size_t i = startPos, count = 0; i < str.size() && count < length; i=i+2) {
        gappedSubstring += str[i];
        ++count;
    }

    return gappedSubstring;
}

std::uint64_t gOMH::gomh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed) {
    if(read.size() < 2*k - 1) return {};
    std::vector<std::uint64_t> hash_vec;
    std::unordered_map<std::string, unsigned> occurrences;
    std::uint64_t cur_seed = seed;

    auto seql = read | seqan3::views::to_char;
    string read_str(seql.begin(), seql.end());

    size_t ll = read_str.size();

    // for(size_t i = 0; i < read_str.size() - k + 1; ++i) {
    //////////////////
    // for(size_t i = 0; i < ll; i += k) {
    //     string kmer;
    //     size_t remaining_length = ll - i;

    //     if (remaining_length >= k) {
    //         kmer = read_str.substr(i, k);
    //     } else if (remaining_length >= k/2){
    //         kmer = read_str.substr(i);
    //     } 
    //     occurrences[kmer]++;
    //     boost::hash_combine(cur_seed, kmer);
    //     boost::hash_combine(cur_seed, occurrences[kmer]);
    //     hash_vec.emplace_back(cur_seed);
    //     cur_seed = seed;
    // }

    for(size_t i = 0; i <= ll - 2*k + 1; ++i) {
        string kmer = getGappedSubstring(read_str, i, k);
        occurrences[kmer]++;
        boost::hash_combine(cur_seed, kmer);
        boost::hash_combine(cur_seed, occurrences[kmer]);
        hash_vec.emplace_back(cur_seed);
        cur_seed = seed;
    }

    auto min_hash = std::min_element(hash_vec.begin(), hash_vec.end());    
    return *min_hash;
}

// using each of all the kmers for bucketing
std::vector<std::uint64_t> gOMH::gomh_pos2(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed) {
    if(read.size() < 2*k - 1) return {};
    std::vector<std::uint64_t> hash_vec;
    std::unordered_map<std::string, unsigned> occurrences;
    std::uint64_t cur_seed = seed;

    auto seql = read | seqan3::views::to_char;
    string read_str(seql.begin(), seql.end());

    size_t ll = read_str.size();

    for(size_t i = 0; i <= ll - 2*k + 1; ++i) {
        string kmer = getGappedSubstring(read_str, i, k);
        occurrences[kmer]++;
        boost::hash_combine(cur_seed, kmer);
        boost::hash_combine(cur_seed, occurrences[kmer]);
        hash_vec.emplace_back(cur_seed);
        cur_seed = seed;
    }
    return hash_vec;
}