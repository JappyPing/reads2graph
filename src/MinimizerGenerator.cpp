// MinimizerGenerator.cpp
#include "MinimizerGenerator.hpp"
#include "ReadWrite.hpp"
#include <algorithm>
#include <execution>
#include <omp.h>
#include <ranges>
#include <boost/functional/hash.hpp>
#include <cmath>

MinimizerGenerator::MinimizerGenerator(cmd_arguments args) : args(args) {}
// MinimizerGenerator::MinimizerGenerator(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : read2count(read2count), args(args) {}

std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator::minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads,std::tuple<unsigned, unsigned, unsigned, double> betterParams)
{   
    // auto better_n = std::get<0>(betterParams);
    auto better_ww = std::get<1>(betterParams);
    auto better_kk = std::get<2>(betterParams);
    // auto prob = std::get<3>(betterParams);
    // auto best_w = round(static_cast<double>(args.read_length) / best_n);   
    auto better_k = static_cast<uint8_t>(better_kk);
    auto better_w = static_cast<uint8_t>(better_ww);

    // int available_cores = omp_get_max_threads();
    // auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    // omp_set_num_threads(num_cores_to_use);

    // #pragma omp parallel for
    // for (size_t i = 0; i < read2count.size(); ++i) {
    //     auto it = std::next(read2count.begin(), i);
    //     const auto& [read, count] = *it;
    // #pragma omp parallel for
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto const & read : unique_reads){
        // auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);
        auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{better_k}) | seqan3::views::minimiser(better_w - better_k + 1);   

        // // Iterate over minimisers and group reads
        for (auto const &minimiser : minimisers) {
            std::uint64_t converted_minimiser = static_cast<std::uint64_t>(minimiser);
            #pragma omp critical
            {
                minimiser_to_reads[converted_minimiser].push_back(read);
            }
        }
    }

    Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Size of minimiser_to_reads: {}.", minimiser_to_reads.size()));  
    return minimiser_to_reads;     
}

// long double MinimizerGenerator::prob(int l, int n, int k, int dt) {
//     int w = round(static_cast<double>(l) / n);
//     // long double p1 = ((w - (round(static_cast<double>(dt) / n) + 1) * k) + 1) / (w - k + 1);
//     long double p1 = (static_cast<long double>(dt) * k) / (n * (w - k + 1));
//     long double p2 = 1 - std::pow(p1, n);
//     return p2;
// }

int MinimizerGenerator::kSize(int L, double p) {
    return ceil((p*(1+L))/(1+p));
}

double MinimizerGenerator::proba(unsigned L, unsigned k) {
    double p;
    p = (static_cast<double>(k))/(L-k+1);
    return p;
}

std::tuple<unsigned, unsigned, unsigned, double> MinimizerGenerator::possibleBetterParameters() {
    unsigned betterK;
    unsigned betterN;
    unsigned betterW;
    double p=0;
    // if (args.read_length >= 8 && args.read_length < 16){
    //     betterK = 3;
    //     betterN = 2;
    //     // betterW = args.read_length;
    //     betterW = round(args.read_length/betterN);
    // } else if (args.read_length >= 16 && args.read_length < 50){
    //     betterK = 4;
    //     betterN = 2;
    //     betterW = round(args.read_length/betterN);
    // } else if (args.read_length >= 50 && args.read_length <= 300) {
    //     // if (args.max_edit_dis == 1 || args.max_edit_dis == 2){
    //     //     betterN = 3;
    //     // } else {
    //     //     betterN = ceil((static_cast<double>(args.max_edit_dis))/2)+1;
    //     // }
    //     betterN = 3;
    //     betterW = round(args.read_length/betterN);
    //     betterK = kSize(betterW, args.bad_kmer_ratio);
    // } 
    if (args.read_length >= 6 && args.read_length < 10){
        betterK = args.k_size;
        betterN = args.window_number;
        betterW = args.read_length;
    } else if (args.read_length >= 10 && args.read_length < 16){
        betterK = args.k_size;
        betterN = args.window_number;
        // betterW = args.read_length;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 16 && args.read_length < 50){
        betterK = args.k_size;
        betterN = args.window_number;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 50 && args.read_length <= 300) {
        // if (args.max_edit_dis == 1 || args.max_edit_dis == 2){
        //     betterN = 3;
        // } else {
        //     betterN = ceil((static_cast<double>(args.max_edit_dis))/2)+1;
        // }
        betterN = args.window_number;
        betterW = round(args.read_length/betterN);
        betterK = kSize(betterW, args.bad_kmer_ratio);
        if (betterK < 4){
            Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Estimated k={} has been changed to 4.", betterK));
            betterK = 4;
        } else 
        if (betterK >= 28) {
            Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Estimated k={} has been changed to 27 as the maximum size of unggaped shape is stricted by 28 in Seqan3.", betterK));  
            betterK = 27;             
        }
    } 
    // if (betterK < 4){
    //     Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Better k {} has been changed to 4.", betterK));
    //     betterK = 4;
    // } else 
    // if (betterK >= 28) {
    //     Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Better k {} has been changed to 27 as the maximum size of unggaped shape is stricted by 28 in Seqan3.", betterK));  
    //     betterK = 27;             
    // }
    if (args.read_length >= 50 && args.read_length <= 300){
        p = 1 - std::pow(proba(betterW, betterK), betterN);
        Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Estimated number of windows: {}, Estimated window size: {}, Estimated K: {} and the probability: {}.", betterN, betterW, betterK, p)); 
    } else {
        Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Number of windows: {}, Window size: {}, K size: {}.", betterN, betterW, betterK));         
    }
    return std::make_tuple(betterN, betterW, betterK, p);   
}

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
//     uint64_t hash_sum = 0;
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

    // auto bestParams = findBestParameters(args.read_length, args.max_edit_dis, args.bad_kmer_ratio);
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

