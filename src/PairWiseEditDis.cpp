#include "PairWiseEditDis.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

PairWiseEditDis::PairWiseEditDis(std::set<std::vector<seqan3::dna5>> const & reads, cmd_arguments args) : reads_(std::move(reads)), args(args) {}

std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> PairWiseEditDis::compute_pairwise_edit_distance()
{
    int min_s = -1 * args.max_edit_dis;
    int max_s = -1 * args.min_edit_dis;

    // Convert the set to a vector
    std::vector<std::vector<seqan3::dna5>> reads_vec(reads_.begin(), reads_.end());

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
        auto it = edge_lst.find(read_pair);
        if (it == edge_lst.end())
        {
            edge_lst[read_pair] = edit_distance;
            edit_distance_counts_[edit_distance]++;            
        }

    }
    std::cout << "Number of total read edges: " << edge_lst.size() << std::endl; 
    for (const auto & [distance, count] : edit_distance_counts_)
    {
        std::cout << "Edit distance by Brute_Force" << distance << ": " << count << " pairs" << std::endl;
    }
    auto filename = args.output_dir / "Brute_Force_pairwise_edit_distance_counts.txt";
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
