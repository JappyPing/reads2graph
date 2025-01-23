// MinimizerGenerator.cpp
#include "MinimizerGenerator.hpp"
#include "ReadWrite.hpp"
#include "omh.hpp"
#include "miniception.hpp"

#include <algorithm>
#include <execution>
#include <omp.h>
#include <ranges>
#include <boost/functional/hash.hpp>
#include <random>
#include <string>
#include <cmath>
#include <boost/format.hpp>
#include <xxhash.hpp>
#include <seeded_prg.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/core/debug_stream.hpp>

MinimizerGenerator::MinimizerGenerator(cmd_arguments args) : args(args) {}

std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator::minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads)
{   
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimiser2reads;
    // #pragma omp parallel for
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto const & read : unique_reads){
        std::vector<std::uint64_t> minimisers;
        uint8_t k_size;
        uint8_t w_size;
        uint8_t num_substr;

        if (args.default_params) {
            if (args.read_length >= 6 && args.read_length < 16){
                args.segmentation = false;
                num_substr = 1;
                k_size = k_estimate(num_substr);  
                w_size = (args.read_length) / 2 - k_size + 1;
                if (w_size == 1){
                    w_size = 2;
                }      
                // w_size = k_size + 1;                   
            } else if (args.read_length >= 16 && args.read_length < 50){
                num_substr = args.substr_number - 1;
            } else if (args.read_length >= 50 && args.read_length <= 300) {
                num_substr = args.substr_number;
            }         
        } else { 
            k_size = args.k_size;
            w_size = args.w_size;
            num_substr = args.substr_number;
        }     
        if (args.read_length >= 6 && args.read_length < 16){ 
            if (args.bucketing_mode == "miniception_gomh") {
                minimisers = Miniception(args).miniception_main(read, k_size, w_size, args.seed);                 
            } else {
                auto minimiser_range = read | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k_size}}) | seqan3::views::minimiser(w_size);    
                std::ranges::copy(minimiser_range, std::back_inserter(minimisers));        
            }
            for (auto const &minimiser : minimisers) {
                #pragma omp critical
                {
                    minimiser2reads[minimiser].push_back(read);
                }
            }                 
        } else {           
            if (args.segmentation){
                auto sub_strs = divide_into_substrings(read, num_substr);
                for (auto const & sub_str : sub_strs){
                    auto substr_size = static_cast<uint8_t>(sub_str.size());
                    if (args.default_params) {
                        k_size = k_estimate(num_substr);
                        if (args.bucketing_mode == "miniception_gomh") {
                            w_size = k_size + 1;
                            // w_size = substr_size - num_substr; // this does not work for miniception
                        } else if (args.bucketing_mode == "minimizer_gomh") {
                            // w_size = k_size + 1;
                            w_size = static_cast<uint8_t>(substr_size * 0.5);
                            // w_size = substr_size - num_substr;
                            // w_size = wSize(k_size, substr_size); 
                        }
                    } 
                    // std::cout << "subread size: " << static_cast<int>(substr_size) << ", k: " << static_cast<int>(k_size) << endl;
                    std::vector<std::uint64_t> minimisers;
                    if (args.bucketing_mode == "miniception_gomh") {
                        minimisers = Miniception(args).miniception_main(sub_str, k_size, w_size, args.seed);                 
                    } else if (args.bucketing_mode == "minimizer_gomh") {
                        auto minimiser_range = sub_str | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k_size}}) | seqan3::views::minimiser(w_size - k_size + 1);              
                        std::ranges::copy(minimiser_range, std::back_inserter(minimisers));        
                    }
                    for (auto const &minimiser : minimisers) {
                        #pragma omp critical
                        {
                            minimiser2reads[minimiser].push_back(read);
                        }
                    }                      
                }
            } else {
                k_size = k_estimate(num_substr);
                uint8_t section_size = static_cast<uint8_t>(std::ceil(args.read_length / num_substr));
                if (args.bucketing_mode == "miniception_gomh") {
                    minimisers = Miniception(args).miniception_main(read, k_size, k_size + 1, args.seed);                 
                } else if (args.bucketing_mode == "minimizer_gomh") {
                    auto minimiser_range = read | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k_size}}) | seqan3::views::minimiser(section_size - k_size + 1);    
                    std::ranges::copy(minimiser_range, std::back_inserter(minimisers));        
                } 
                for (auto const &minimiser : minimisers) {
                    #pragma omp critical
                    {
                        minimiser2reads[minimiser].push_back(read);
                    }
                }                        
            }
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, boost::str(boost::format("Size of minimiser2reads: %1%!") % minimiser2reads.size()));
    return minimiser2reads;     
}

// Function to split each sequence in the read vector into substrings
std::vector<std::vector<seqan3::dna5>> MinimizerGenerator::divide_into_substrings(const std::vector<seqan3::dna5> & read, int num_substrs)
{
    std::vector<std::vector<seqan3::dna5>> substrings;
    int total_length = read.size();
    int window_size = total_length / num_substrs;

    for (int i = 0; i < num_substrs; ++i)
    {
        int start = i * window_size;
        int end = (i == num_substrs - 1) ? total_length : start + window_size;

        // Slice the sequence for this window and store it in a vector
        auto subrange = read | seqan3::views::slice(start, end);
        substrings.emplace_back(subrange.begin(), subrange.end());
    }

    return substrings;
}

// int MinimizerGenerator::kSize(int L, double p) {
//     return ceil((p*(1+L))/(1+p));
// }

uint8_t MinimizerGenerator::wSize(uint8_t k, uint8_t read_len) {
    uint8_t w;
    if (args.bucketing_mode == "miniception_gomh") {
        w = k + 1;
    } else if (args.bucketing_mode == "minimizer_gomh") {
        // w = static_cast<uint8_t>(std::ceil(std::pow(4, k / 4)));
        if ((w > read_len) || (w < k)){
            w = 2*k;
        }
    }
    return w;
}

// uint8_t MinimizerGenerator::k_estimate(uint8_t num_substr, uint8_t read_size) {
//     double L_seg = static_cast<double>(read_size) / num_substr;
//     double term = std::pow(1 - args.probability, 1.0 / num_substr);
//     double k_max = std::min((term * (L_seg + 1)) / (args.max_edit_dis / static_cast<double>(num_substr) + term), L_seg / args.max_edit_dis);
//     uint8_t k = std::min(std::max(4, static_cast<int>(std::floor(k_max))), 28);
//     return k;
// }

uint8_t MinimizerGenerator::k_estimate(uint8_t N) {
    int segment_size = args.read_length / N;
    uint8_t k = static_cast<uint8_t>(ceil((args.differ_kmer_ratio * N * (1 + segment_size))/(2 + N * args.differ_kmer_ratio)));
    if (k > 28) {
        k = 28;       
    } else if (k < 4){
        k = 4;       
    }
    return k;
}

// uint8_t MinimizerGenerator::k_estimate(uint8_t N) {
//     uint8_t k;
//     if (args.segmentation) {
//         int segment_size = args.read_length / N;
//         k = static_cast<uint8_t>(ceil((args.differ_kmer_ratio * N * (1 + segment_size))/(2 + N * args.differ_kmer_ratio)));
//     } else {
//         // p = (L-k+1 - d*k)/(L-k+1)
//         k = ceil((1-args.probability)*(1+args.read_length)/(1+args.max_edit_dis-args.probability));         
//     }
//     if (k > 28) {
//         k = 28;       
//     } else if (k < 4){
//         k = 4;       
//     }
//     return k;
// }


double MinimizerGenerator::proba(unsigned L, unsigned k) {
    double p;
    p = (static_cast<double>(k))/(L-k+1);
    return p;
}

/*
std::tuple<unsigned, unsigned, unsigned, double> MinimizerGenerator::k_w_estimate() {
    unsigned betterK;
    unsigned betterN;
    unsigned betterW;
    double p=0;
    if (args.read_length >= 6 && args.read_length < 16){
        betterK = args.k_size;
        betterN = 1;
        betterW = args.read_length;
    } else if (args.read_length >= 16 && args.read_length < 50){
        betterK = args.k_size;
        betterN = 2;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 50 && args.read_length <= 300) {
        betterN = args.substr_number;
        betterW = round(args.read_length/betterN);
        betterK = kSize(betterW, args.differ_kmer_ratio);
        if (betterK < 4){
            // Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Estimated k={} has been changed to 4.", betterK));
            Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("Estimated k=%1% has been changed to 4.") % betterK));

            betterK = 4;
        } else if (betterK >= 28) {
            Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("Estimated k=%1% has been changed to 27 as the maximum size of unggaped shape is stricted by 28 in Seqan3.") % betterK)); 
            betterK = 27;             
        }
    } 
    if (args.read_length >= 50 && args.read_length <= 300){
        p = 1 - std::pow(proba(betterW, betterK), betterN);
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("Estimated number of windows: %1%, Estimated window size: %2%, Estimated K: %3% and the probability: %4%.") % betterN % betterW % betterK % p)); 
    } else {
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("Number of windows: %1%, Window size: %2%, K size: %3%.") % betterN % betterW % betterK));         
    }
    auto number_kmer = betterW - betterK + 1;
    if ((number_kmer) < 3 ){
        Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("only %1% kmers setted for minimizer selection.") % number_kmer));
    }

    return std::make_tuple(betterN, betterW, betterK, p);   
}

std::tuple<unsigned, unsigned, unsigned, double> MinimizerGenerator::overlappingWindowParameters() {
    unsigned betterK;
    unsigned betterN;
    unsigned betterW;
    double p=0;
    if (args.read_length >= 6 && args.read_length < 10){
        betterK = args.k_size;
        betterN = args.substr_number;
        betterW = args.read_length;
    } else if (args.read_length >= 10 && args.read_length < 16){
        betterK = args.k_size;
        betterN = args.substr_number;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 16 && args.read_length < 50){
        betterK = args.k_size;
        betterN = args.substr_number;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 50 && args.read_length <= 300) {
        betterN = args.substr_number;
        betterW = round(args.read_length/betterN);
        betterK = kSize(betterW, args.differ_kmer_ratio);
        if (betterK < 4){
            Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("Estimated k=%1% has been changed to 4.") % betterK));

            betterK = 4;
        } else if (betterK >= 28) {
            Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("Estimated k=%1% has been changed to 27 as the maximum size of unggaped shape is stricted by 28 in Seqan3.") % betterK)); 
            betterK = 27;             
        }
    } 
    if (args.read_length >= 50 && args.read_length <= 300){
        p = 1 - std::pow(proba(betterW, betterK), betterN);
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("Estimated number of windows: %1%, Estimated window size: %2%, Estimated K: %3% and the probability: %4%.") % betterN % betterW % betterK % p)); 
    } else {
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("Number of windows: %1%, Window size: %2%, K size: %3%.") % betterN % betterW % betterK));         
    }
    auto number_kmer = betterW - betterK + 1;
    if ((number_kmer) < 3 ){
        Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("only %1% kmers setted for minimizer selection.") % number_kmer));
    }

    return std::make_tuple(betterN, betterW, betterK, p);   
}

// using minimizer only with random ordering to bucket reads multiple times
std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator::minimizer_only2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, uint8_t k, unsigned m)
{
    std::mt19937_64 generator(args.seed);
    // Specify the range of values for your seeds
    std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());
    for (unsigned i = 0; i < m; ++i) {
        // Generate a new seed for this hash iteration
        std::uint64_t cur_seed = distribution(generator);
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (auto const & read : unique_reads){
            unsigned int len_read = read.size();
            // Compute the minimizer hash for the entire read (produces one minimizer)
            auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{k}}, seqan3::window_size{len_read - k + 1}, seqan3::seed{cur_seed});

            // Get the first (and only) minimizer since the window spans the entire read
            std::uint64_t converted_minimiser = static_cast<std::uint64_t>(*minimisers.begin());
            #pragma omp critical
            {
                minimiser_to_reads[converted_minimiser].push_back(read);
            }
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, boost::str(boost::format("Size of minimiser_to_reads: %1%!") % minimiser_to_reads.size()));
    return minimiser_to_reads;     
}
*/

// // Function to sample p percentage of elements
// std::size_t MinimizerGenerator::min_k_in_sampling(const std::vector<std::vector<seqan3::dna5>>& unique_reads, double p) {
//     auto total_uniq_num = unique_reads.size();
//     Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of unique reads: {} ", total_uniq_num));
//     std::size_t sample_num = static_cast<std::size_t>(p * total_uniq_num);

//     // Set up a default random number generator
//     std::default_random_engine gen;

//     // Create an index vector to randomly shuffle
//     std::vector<std::size_t> indices(total_uniq_num);
//     std::iota(indices.begin(), indices.end(), 0);

//     // Shuffle the indices
//     std::shuffle(indices.begin(), indices.end(), gen);

//     // Create a vector to store the sampled sequences
//     std::vector<std::vector<seqan3::dna5>> sampled_uniq_reads;

//     std::size_t min_length = std::numeric_limits<std::size_t>::max();
//     // Extract the first num_to_sample indices from the shuffled vector
//     #pragma omp parallel for
//     for (std::size_t i = 0; i < sample_num; ++i) {
//         #pragma omp critical
//         {
//            sampled_uniq_reads.push_back(unique_reads[indices[i]]); 
//            auto cur_read_length = unique_reads[indices[i]].size();
//             if (cur_read_length < min_length) {
//                 min_length = cur_read_length;
//             }
//         }        
//     }
//     auto max_k = round(min_length/3)-1;
//     // Parallelize the loop using OpenMP
//     // #pragma omp parallel for
//     for (std::size_t k = 3; k < max_k; ++k) {
//         bool all_unique_flag = true;
//         for (auto& cur_seq : sampled_uniq_reads) {
//             auto kmer_view = cur_seq | seqan3::views::kmer_hash(seqan3::ungapped{k});
//             std::set<size_t> kmer_set(kmer_view.begin(), kmer_view.end());

//             if (kmer_set.size() < (cur_seq.size() - k + 1)){
//                 all_unique_flag = false;
//                 break;
//             }
//         }
//         if (all_unique_flag){
//             return k;
//         }
//     }
//     return max_k - 1;
// }

// splicing multiple minimisers as one key
// void MinimizerGenerator::process_read(const std::vector<seqan3::dna5> &read)
// {
//     // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
//     // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
//     auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);
//     std::uint64_t hash_sum = 0;
//     for (auto const &minimiser : minimisers) {
//         // std::uint64_t converted_minimiser = static_cast<std::uint64_t>(minimiser);
//         // seqan3::debug_stream << "Minimiser: " << minimiser << " " << converted_minimiser << '\n';
//         hash_sum += minimiser;
//     } 
//     minimiser_to_reads[hash_sum].push_back(read);
// }

// std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator::process_reads_in_parallel()
// {   
    // size_t max_length = 0;
    // size_t min_length = std::numeric_limits<size_t>::max(); // Set to maximum possible value initially
    // // for (const auto& read : unique_reads_) {
    // for (const auto& [read, count] : read2count) {
    //     size_t read_length = read.size();
    //     if (read_length > max_length) {
    //         max_length = read_length;
    //     }
    //     // Update minimum length
    //     if (read_length < min_length) {
    //         min_length = read_length;
    //     }
    // }    
    // std::cout << "Minimum read Length: " << min_length << endl;
    // int best_k;
    // int best_w;

    // auto bestParams = findBestParameters(args.read_length, args.max_edit_dis, args.differ_kmer_ratio);
    // auto best_n = std::get<0>(bestParams);
    // auto best_kk = std::get<1>(bestParams);
    // auto best_w = round(static_cast<double>(args.read_length) / best_n);
    // // std::cout << "Best number of windows: " << best_n << "\n" << "Best K: " << best_kk << "\n" << "Best probability: " << std::get<2>(bestParams) << std::endl;   
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Best number of windows: {}, Best K: {} and Best probability: {}.", best_n, best_kk, std::get<2>(bestParams)));     
    // auto best_k = static_cast<uint8_t>(best_kk);

    // int desired_num_cores = 26; /* specify the number of cores you want to use */
    // int available_cores = omp_get_max_threads();
    // Declare and define a global variable for available cores
    // int available_cores = omp_get_max_threads();
    // Ensure the user-specified number of cores is within a valid range
    // int num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);

    // Set the number of threads for OpenMP
    // omp_set_num_threads(num_cores_to_use);
    // auto w_size = args.window_size;
    
    // for (size_t i = 0; i < unique_reads_.size(); ++i)
    // {
    //     auto it = std::next(unique_reads_.begin(), i);
    //     // Dereference the iterator to get the actual read content
    //     auto const &read = *it;

    // int available_cores = omp_get_max_threads();
    // auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    // omp_set_num_threads(num_cores_to_use);

    // #pragma omp parallel for
    // for (size_t i = 0; i < read2count.size(); ++i) {
    //     auto it = std::next(read2count.begin(), i);
    //     const auto& [read, count] = *it;
    //     // auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);
    //     auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{best_k}) | seqan3::views::minimiser(best_w - best_k + 1);   
    
        // forward and backward minimisers
        // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{args.k_size}}, seqan3::window_size{w_size - args.k_size + 1});


        // Iterate over minimisers and group reads
        // std::size_t seed = 0;
        // for (auto const &minimiser : minimisers) {
        //     std::uint64_t converted_minimiser = static_cast<std::uint64_t>(minimiser);
        //     boost::hash_combine(seed, converted_minimiser);
        // }
        // #pragma omp critical
        // {
        //     minimiser_to_reads[seed].push_back(read);
        // }
        // // Iterate over minimisers and group reads
    //     for (auto const &minimiser : minimisers) {
    //         std::uint64_t converted_minimiser = static_cast<std::uint64_t>(minimiser);
    //         #pragma omp critical
    //         {
    //             minimiser_to_reads[converted_minimiser].push_back(read);
    //         }
    //     }
    // }
    // #pragma omp taskwait
    // #pragma omp parallel for
    // for (size_t i = 0; i < unique_reads_.size(); ++i)
    // {
    //     auto it = std::next(unique_reads_.begin(), i);
    //     // Process each read in parallel
    //     #pragma omp task
    //     {
    //         process_read(*it);
    //     }   
    // }
    // #pragma omp taskwait
    // std::for_each(std::execution::par_unseq, unique_reads_.begin(), unique_reads_.end(),
    //                 [this](const auto &read)
    //                 {
    //                     // Process each read in parallel
    //                     process_read(read);
    //                 });

    // // Erase elements where the size of the associated vector is equal to 1
    // for (auto it = minimiser_to_reads.begin(); it != minimiser_to_reads.end(); )
    // {
    //     if (it->second.size() == 1)
    //     {
    //         it = minimiser_to_reads.erase(it);
    //     }
    //     else
    //     {
    //         ++it;
    //     }
    // }

    // for (auto const &[minimiser, reads] : minimiser_to_reads) {
    //     // seqan3::debug_stream << "Minimiser for current Block: " << minimiser << '\n';
    //     seqan3::debug_stream << "Reads in current Block:\n";
    //     // for (auto const &read : reads) {
    //     //     seqan3::debug_stream << read << '\n';
    //     // }
    //     std::cout << "Number of reads: " << reads.size() << std::endl;
    //     // seqan3::debug_stream << "-----------------\n";
    //     // Utils::getInstance().logger(LOG_LEVEL_INFO,  "-----------------");
    //     // break;
    // } 
    // std::cout << "Size of minimiser_to_reads: " << minimiser_to_reads.size() << std::endl;  
//     Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Size of minimiser_to_reads: {}.", minimiser_to_reads.size()));  
//     return minimiser_to_reads;     
// }

// long double MinimizerGenerator::prob(int l, int n, int k, int dt) {
//     int w = round(static_cast<double>(l) / n);
//     // long double p1 = ((w - (round(static_cast<double>(dt) / n) + 1) * k) + 1) / (w - k + 1);
//     long double p1 = (static_cast<long double>(dt) * k) / (n * (w - k + 1));
//     long double p2 = 1 - std::pow(p1, n);
//     return p2;
// }