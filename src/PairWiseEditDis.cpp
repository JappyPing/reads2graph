#include "PairWiseEditDis.hpp"
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

PairWiseEditDis::PairWiseEditDis(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : read2count(read2count), args(args) {}

// std::unordered_map<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>, int, unordered_pair> PairWiseEditDis::compute_pairwise_edit_distance()
std::map<std::set<std::vector<seqan3::dna5>>, int> PairWiseEditDis::compute_pairwise_edit_distance()

{
    // Convert the set to a vector
    // std::vector<std::vector<seqan3::dna5>> reads_vec(reads_.begin(), reads_.end());
    std::vector<std::vector<seqan3::dna5>> reads_vec;
    std::transform(read2count.begin(), read2count.end(),
                   std::back_inserter(reads_vec),
                   [](const auto& pair) { return pair.first; });

    int min_s = -1 * args.max_edit_dis;
    int max_s = -1 * args.min_edit_dis;

    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s} | seqan3::align_cfg::output_score{};

    auto pairwise_combinations = seqan3::views::pairwise_combine(reads_vec);

    // Declare and define a global variable for available cores
    // int available_cores = omp_get_max_threads();
    // // std::cout << "The maximum number of threads available:" << available_cores << std::endl;
    // // Ensure the user-specified number of cores is within a valid range
    // int num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);

    // // Set the number of threads for OpenMP
    // omp_set_num_threads(num_cores_to_use);

    #pragma omp parallel for
    for (size_t i = 0; i < pairwise_combinations.size(); ++i)
    {
        auto const &combination = pairwise_combinations[i];
        auto const &seq1 = std::get<0>(combination);
        auto const &seq2 = std::get<1>(combination);

        std::set<std::vector<seqan3::dna5>> read_pair_set;
        read_pair_set.insert(seq1);
        read_pair_set.insert(seq2);
        auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
        // Iterate over alignment results and access the scores
        for (auto const &result : alignment_results)
        {
            int edit_distance = result.score();
            if ((edit_distance >= min_s) && (edit_distance <= max_s)) 
            {
                #pragma omp critical
                {
                    edge_lst[read_pair_set] = edit_distance;
                    edit_distance_counts_[edit_distance]++; 
                }                    
            }
        }

        /*
        auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
        std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>> read_pair;
        if (seq1 == seq2) {
            throw std::runtime_error("Error: two reads are the same in the unoredered read pair!");
        } else if (seq1 < seq2) {
            auto read_pair = std::make_pair(seq1, seq2);
        } else {
            auto read_pair = std::make_pair(seq2, seq1);       
        }
        // auto it = edge_lst.find(read_pair);
        // Iterate over alignment results and access the scores
        for (auto const &result : alignment_results)
        {
            int edit_distance = result.score();
            if ((edit_distance >= min_s) && (edit_distance <= max_s)) 
            {
                #pragma omp critical
                {
                    edge_lst[read_pair] = edit_distance;
                    edit_distance_counts_[edit_distance]++; 
                }                    
            }
        }
        */
    }        

    // std::cout << "Number of total read edges: " << edge_lst.size() << std::endl; 
    Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Number of total read edges: {}.", edge_lst.size()));
    for (const auto & [distance, count] : edit_distance_counts_)
    {
        // std::cout << "Edit distance by Brute_Force" << distance << ": " << count << " pairs" << std::endl;
        Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Edit distance by Brute_Force {}: {} pairs.", distance, count));
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

/*
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
        edge_lst[read_pair] = edit_distance;
        edit_distance_counts_[edit_distance]++; 

        // auto it = edge_lst.find(read_pair);
        // if (it == edge_lst.end())
        // {
        //     edge_lst[read_pair] = edit_distance;
        //     edit_distance_counts_[edit_distance]++;            
        // }

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
*/
 
