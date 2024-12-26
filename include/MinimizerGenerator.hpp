// MinimizerGenerator.hpp

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <vector>
#include <set>
#include <seqan3/alphabet/all.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits> // for std::decay_t
#include "Utils.hpp"
#include "LoggingLevels.hpp"

using namespace std;
using namespace seqan3::literals;

class MinimizerGenerator
{
public:
    MinimizerGenerator(cmd_arguments args);
    // MinimizerGenerator(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args);
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::tuple<unsigned, unsigned, unsigned, double> betterParams);
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads);
    double proba(unsigned L, unsigned k);
    int kSize(int L, double p);
    std::vector<std::vector<seqan3::dna5>> divide_into_substrings(const std::vector<seqan3::dna5> & read, int num_substrs);
    std::tuple<unsigned, unsigned, unsigned, double> overlappingWindowParameters();
    std::tuple<unsigned, unsigned, unsigned, double> segmentation_parameters();
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer_only2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, uint8_t k, unsigned m);
    // std::tuple<int, int, double> findBestParameters(int l, int dt, double pt);
    // long double prob(int l, int n, int k, int dt);
    // void process_read(const std::vector<seqan3::dna5> &read);
    // std::size_t min_k_in_sampling(const std::vector<std::vector<seqan3::dna5>>& unique_reads, double p);

private:
    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads;
    cmd_arguments args;
};
    // Function to group reads into blocks based on stored minimiser sets
    // void group_reads_into_blocks()
    // {
    //     std::unordered_map<decltype(minimiser_sets_), std::vector<std::vector<seqan3::dna5>>> blocks;

    //     for (const auto &read : unique_reads_)
    //     {
    //         auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(5);

    //         // Find the corresponding minimiser set in the stored sets
    //         auto set_it = minimiser_sets_.find(minimisers);

    //         if (set_it != minimiser_sets_.end())
    //         {
    //             // Add the read to the corresponding block
    //             blocks[*set_it].push_back(read);
    //         }
    //     }

    //     // Process the grouped blocks as needed
    //     for (const auto &[minimiser_set, reads_in_block] : blocks)
    //     {
    //         std::cout << "Minimiser Set: " << minimiser_set << '\n';
    //         std::cout << "Reads in Block:\n";
    //         for (const auto &read : reads_in_block)
    //         {
    //             std::cout << read << '\n';
    //         }
    //         std::cout << "-----------------\n";
    //     }
    // }

// class MinimizerGenerator
// {
// public:
//     MinimizerGenerator(std::set<std::vector<seqan3::dna5>> unique_reads) : unique_reads_(std::move(unique_reads)) {}

//     void process_read(const std::vector<seqan3::dna5> &read)
//     {
//         // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
//         // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
//         auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(5);

//         seqan3::debug_stream << read << '\n';
//         seqan3::debug_stream << minimisers << '\n';
//     }

//     void process_reads_in_parallel()
//     {
//         std::for_each(std::execution::par_unseq, unique_reads_.begin(), unique_reads_.end(),
//                       [this](const auto &read)
//                       {
//                           // Process each read in parallel
//                           process_read(read);
//                       });
//     }

// private:
//     std::set<std::vector<seqan3::dna5>> unique_reads_;
// };


// class MinimizerGenerator
// {
// public:
//     MinimizerGenerator(std::set<std::vector<seqan3::dna5>> unique_reads) : unique_reads_(std::move(unique_reads)) {}

//     void process_read(const std::vector<seqan3::dna5> & read)
//     {
//         // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
//         // auto minimisers = read | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
//         auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(5);

//         seqan3::debug_stream << read << '\n';
//         seqan3::debug_stream << minimisers << '\n';
//     }

//     void process_reads_in_parallel()
//     {
//         // Enable parallel execution
//         #pragma omp parallel for
//         for (auto it = unique_reads_.begin(); it != unique_reads_.end(); ++it)
//         {
//             // Process each read in parallel
//             process_read(*it);
//         }
//     }

// private:
//     std::set<std::vector<seqan3::dna5>> unique_reads_;
// };