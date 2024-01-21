// MinimizerGenerator.cpp
#include "MinimizerGenerator.h"
#include <algorithm>
#include <execution>
#include <omp.h>
#include <ranges>

MinimizerGenerator::MinimizerGenerator(std::set<std::vector<seqan3::dna5>> unique_reads, cmd_arguments args) : unique_reads_(std::move(unique_reads)), args(args) {}

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
    for (auto const &minimiser : minimisers) {
        std::int64_t converted_minimiser = static_cast<std::int64_t>(minimiser);
        // seqan3::debug_stream << "Minimiser: " << minimiser << " " << converted_minimiser << '\n';
        minimiser_to_reads[converted_minimiser].push_back(read);
    }
}

std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator::process_reads_in_parallel()
{
    // int desired_num_cores = 26; /* specify the number of cores you want to use */
    // int available_cores = omp_get_max_threads();
    // Declare and define a global variable for available cores
    int available_cores = omp_get_max_threads();
    // Ensure the user-specified number of cores is within a valid range
    int num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);

    // Set the number of threads for OpenMP
    omp_set_num_threads(num_cores_to_use);

    #pragma omp parallel for
    for (size_t i = 0; i < unique_reads_.size(); ++i)
    {
        auto it = std::next(unique_reads_.begin(), i);
        // Dereference the iterator to get the actual read content
        auto const &read = *it;

        auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) |
                        seqan3::views::minimiser(args.window_size - args.k_size + 1);

        // Iterate over minimisers and group reads
        for (auto const &minimiser : minimisers) {
            std::int64_t converted_minimiser = static_cast<std::int64_t>(minimiser);
            #pragma omp critical
            {
                minimiser_to_reads[converted_minimiser].push_back(read);
            }
        }
    }
    #pragma omp taskwait
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

    // Erase elements where the size of the associated vector is equal to 1
    for (auto it = minimiser_to_reads.begin(); it != minimiser_to_reads.end(); )
    {
        if (it->second.size() == 1)
        {
            it = minimiser_to_reads.erase(it);
        }
        else
        {
            ++it;
        }
    }

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