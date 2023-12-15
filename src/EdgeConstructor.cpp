// EdgeConstructor.cpp
#include "EdgeConstructor.h"
#include "GraphManager.h"
#include <algorithm>
#include <ranges>
#include <utility>
#include <vector>
#include <execution>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

EdgeConstructor::EdgeConstructor(std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads, cmd_arguments args) : minimiser_to_reads_(std::move(minimiser_to_reads)), args(args) {}

void EdgeConstructor::process_block(const std::vector<std::vector<seqan3::dna5>> &reads_vec, int min_s, int max_s)
{
    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s}
                    | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_sequence1_id{}
                    | seqan3::align_cfg::output_sequence2_id{};

        auto alignment_results = seqan3::align_pairwise(seqan3::views::pairwise_combine(reads_vec), config);
        auto filter_v = std::views::filter(
            [&min_s, &max_s](auto && res)
            {
                return (res.score() > min_s) && (res.score() <= max_s);
            });

    // Create a mapping, that maps from alignment-pair-ids → sequence-ids
    auto pair_indices = seqan3::views::pairwise_combine(std::views::iota(0ul, reads_vec.size()));
    for (auto const & result : alignment_results | filter_v)  {
        // seqan3::debug_stream << result << '\n';

        // The sequence1_id is reporting the index of the alignment-pair, not the ids of the pairs them self
        auto alignment_id = result.sequence1_id(); // result.sequence1_id() == result.sequence2_id()

        auto [sid1, sid2] = pair_indices[alignment_id]; // alignment-pair-id -> sequence-id
        // seqan3::debug_stream << "sid1: " << sid1 << "; sid2: " << sid2 << "\n"; // the real sequence ids
        auto read_pair = std::make_pair(reads_vec[sid1], reads_vec[sid2]);
        int edit_distance = result.score();
        edge_lst[read_pair] = edit_distance;
    }
}

void EdgeConstructor::process_blocks_in_parallel()
{
    int min_s = -1 * args.max_edit_dis;
    int max_s = -1 * args.min_edit_dis;
    std::for_each(std::execution::par_unseq, minimiser_to_reads_.begin(), minimiser_to_reads_.end(),
                    [this, min_s, max_s](const auto &entry) {
                        // std::int64_t minimiser = entry.first;
                        const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;
                        // Process each group in parallel
                        process_block(reads_vec, min_s, max_s);
                    });
}

std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> EdgeConstructor::get_edge_lst() 
{
    process_blocks_in_parallel();
    // Print the stored read pairs and edit distances
    for (const auto &[read_pair, edit_distance] : edge_lst)
    {
        edit_distance_counts_[edit_distance]++;
        // seqan3::debug_stream << "Read Pair: " << read_pair.first << " - " << read_pair.second
                            // << ", Edit Distance: " << edit_distance << '\n';
    }   
    std::cout << "Number of read edges: " << edge_lst.size() << std::endl;
    
    for (const auto & [distance, count] : edit_distance_counts_)
    {
        std::cout << "Edit distance by minimiser" << distance << ": " << count << " pairs" << std::endl;
    }
    auto filename = args.output_dir / "Minimiser_edit_distance_counts.txt";
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (const auto & [distance, count] : edit_distance_counts_)
    {
        file << "Edit distance " << distance << ": " << count << " pairs\n";
    }

    file.close();

    return edge_lst;
}

// This method get the read pair firstly, and then calculate the read pair edit distance one by one in parallel. It runs slowly.
// void EdgeConstructor::process_block()
// {
//     // int maxInt = std::numeric_limits<int>::max();
//     // Configure the alignment kernel.
    
//     int min_s = -1 * args.max_edit_dis;
//     for (const auto &entry : minimiser_to_reads_) 
//     {
//         // Create a map to store read pairs and edit distances
//         std::map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, pair_comparator> edge_lst;

//         // Function to calculate edit distance and store in the map
//         auto calculate_and_store_edit_distance = [&](auto && pair)
//         {
//             auto &seq1 = std::get<0>(pair);
//             auto &seq2 = std::get<1>(pair);

//             auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

//             auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);

//             for (auto &&result : alignment_results)
//             {
//                 int edit_distance = result.score();
//                 if (edit_distance >= min_s)
//                 {
//                     auto read_pair = std::make_pair(seq1, seq2);
//                     auto it = edge_lst.find(read_pair);

//                     // Check if the pair is not in the map or has a smaller edit distance
//                     if (it == edge_lst.end() || edit_distance < it->second)
//                     {
//                         edge_lst[read_pair] = edit_distance;
//                     }
//                 }
//             }
//         };

//         // std::int64_t minimiser = entry.first;
//         const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;
//         // Parallelize the calculation using C++20 parallel algorithms
//         auto pair_vec = seqan3::views::pairwise_combine(reads_vec);
//         std::for_each(std::execution::par_unseq, pair_vec.begin(), pair_vec.end(), calculate_and_store_edit_distance);

//         std::cout << "Number of read edges: " << edge_lst.size() << std::endl; 
//         // Print the stored read pairs and edit distances
//         // for (const auto &[read_pair, edit_distance] : edge_lst)
//         // {
//         //     seqan3::debug_stream << "Read Pair: " << read_pair.first << " - " << read_pair.second
//         //                         << ", Edit Distance: " << edit_distance << '\n';
//         // }
//     }

 
// }