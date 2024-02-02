// MinimizerGenerator.cpp
#include "MinimizerGenerator.h"
#include "GraphManager.h"
#include "ReadWrite.h"
#include <algorithm>
#include <execution>
#include <omp.h>
#include <ranges>
#include <boost/functional/hash.hpp>
#include <cmath>

// MinimizerGenerator::MinimizerGenerator(std::set<std::vector<seqan3::dna5>> unique_reads, cmd_arguments args) : unique_reads_(std::move(unique_reads)), args(args) {}
MinimizerGenerator::MinimizerGenerator(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : read2count(read2count), args(args) {}

// splicing multiple minimisers as one key
// void MinimizerGenerator::process_read(const std::vector<seqan3::dna5> &read)
// {
//     // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
//     // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
//     auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);
//     uint64_t hash_sum = 0;
//     for (auto const &minimiser : minimisers) {
//         // std::int64_t converted_minimiser = static_cast<std::int64_t>(minimiser);
//         // seqan3::debug_stream << "Minimiser: " << minimiser << " " << converted_minimiser << '\n';
//         hash_sum += minimiser;
//     } 
//     minimiser_to_reads[hash_sum].push_back(read);
// }
/*
void MinimizerGenerator::process_read(const std::vector<seqan3::dna5> &read)
{
    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
    // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
    auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);

    // seqan3::debug_stream << read << '\n';
    // seqan3::debug_stream << minimisers << '\n';
    // Iterate over the elements of the minimiser range
    // for (auto const &minimiser : minimisers)
    // {
    //     // Process each minimiser as needed
    //     seqan3::debug_stream << "Minimiser: " << minimiser << '\n';
    // }

    // Iterate over minimisers and group reads
    // using each minimiser as a key
    for (auto const &minimiser : minimisers) {
        std::int64_t converted_minimiser = static_cast<std::int64_t>(minimiser);
        // seqan3::debug_stream << "Minimiser: " << minimiser << " " << converted_minimiser << '\n';
        minimiser_to_reads[converted_minimiser].push_back(read);
    } 
}
*/

std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator::process_reads_in_parallel()
{
    size_t max_length = 0;
    size_t min_length = std::numeric_limits<size_t>::max(); // Set to maximum possible value initially
    // for (const auto& read : unique_reads_) {
    for (const auto& [read, count] : read2count) {
        size_t read_length = read.size();
        if (read_length > max_length) {
            max_length = read_length;
        }
        // Update minimum length
        if (read_length < min_length) {
            min_length = read_length;
        }
    }    
    std::cout << "Minimum read Length: " << min_length << endl;
    auto bestParams = findBestParameters(min_length, args.max_edit_dis);
    auto best_n = std::get<0>(bestParams);
    auto best_k = std::get<1>(bestParams);
    auto best_w = static_cast<uint8_t>(round(static_cast<double>(min_length) / best_n));
    std::cout << "Best number of windows: " << best_n << "\n" << "Best K: " << best_k << "\n" << "Best probability: " << std::get<2>(bestParams) << std::endl;

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
    #pragma omp parallel for
    for (size_t i = 0; i < read2count.size(); ++i) {
        auto it = std::next(read2count.begin(), i);
        const auto& [read, count] = *it;
        // auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);
        auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{best_k}) | seqan3::views::minimiser(best_w - best_k + 1);        
        // forward and backward minimisers
        // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{args.k_size}}, seqan3::window_size{w_size - args.k_size + 1});


        // Iterate over minimisers and group reads
        // std::size_t seed = 0;
        // for (auto const &minimiser : minimisers) {
        //     std::int64_t converted_minimiser = static_cast<std::int64_t>(minimiser);
        //     boost::hash_combine(seed, converted_minimiser);
        // }
        // #pragma omp critical
        // {
        //     minimiser_to_reads[seed].push_back(read);
        // }
        // // Iterate over minimisers and group reads
        for (auto const &minimiser : minimisers) {
            std::int64_t converted_minimiser = static_cast<std::int64_t>(minimiser);
            #pragma omp critical
            {
                minimiser_to_reads[converted_minimiser].push_back(read);
            }
        }
    }
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
    //     // Utils::logMessage(LOG_LEVEL_INFO,  "-----------------");
    //     // break;
    // } 
    std::cout << "Size of minimiser_to_reads: " << minimiser_to_reads.size() << std::endl;  
    return minimiser_to_reads;     
}

long double MinimizerGenerator::prob(int l, int n, int k, int dt) {
    int w = round(static_cast<double>(l) / n);
    // long double p1 = ((w - (round(static_cast<double>(dt) / n) + 1) * k) + 1) / (w - k + 1);
    long double p1 = (static_cast<long double>(dt) * k) / (n * (w - k + 1));
    long double p2 = 1 - std::pow(p1, n);
    return p2;
}

std::tuple<int, int, double> MinimizerGenerator::findBestParameters(int l, int dt) {
    int bestN = 1;
    int bestK = 3;
    long double bestResult = 0.0;
    
    for (int n = 1; n <= round(dt/2); ++n) {
        // int max_k = round(static_cast<double>(l) / n);
        int max_k = round(l / n);
        for (int k = 3; k < max_k; ++k) {
            auto k4 = std::pow(4, k);
            if (k4 >= (l - k + 1) * 10000){
                long double result = prob(l, n, k, dt);
                if (result > bestResult) {
                    bestResult = result;
                    bestN = n;
                    bestK = k;
                }                
            } else{
                continue;
            }

        }
    }

    return std::make_tuple(bestN, bestK, bestResult);
}