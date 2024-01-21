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
#include <seqan3/alphabet/views/to_rank.hpp>
#include <omp.h>
#include <eigen3/Eigen/Dense>
// #include <thread>
// #include <mutex>

EdgeConstructor::EdgeConstructor(std::unordered_map<std::int64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads, cmd_arguments args) : minimiser_to_reads_(std::move(minimiser_to_reads)), args(args) {}

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

void EdgeConstructor::process_blocks_in_parallel()
{
    int min_s = -1 * args.max_edit_dis;
    int max_s = -1 * args.min_edit_dis;

    // Declare and define a global variable for available cores
    int available_cores = omp_get_max_threads();
    // Ensure the user-specified number of cores is within a valid range
    int num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);

    // Set the number of threads for OpenMP
    omp_set_num_threads(num_cores_to_use);

    #pragma omp parallel for
    for (auto i = 0u; i < minimiser_to_reads_.size(); ++i) {
        const auto &entry = *std::next(minimiser_to_reads_.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;

        // Use OpenMP to parallelize the processing of each group of reads
        #pragma omp task
        {
            // process_block(reads_vec, min_s, max_s);
            auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s} | seqan3::align_cfg::output_score{};

            auto pairwise_combinations = seqan3::views::pairwise_combine(reads_vec);

            #pragma omp parallel for
            for (size_t i = 0; i < pairwise_combinations.size(); ++i)
            {
                auto const &combination = pairwise_combinations[i];
                auto const &seq1 = std::get<0>(combination);
                auto const &seq2 = std::get<1>(combination);

                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);

                auto read_pair = std::make_pair(seq1, seq2);
                auto it = edge_lst.find(read_pair);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results)
                {
                    int edit_distance = result.score();
                    if (it == edge_lst.end() && (edit_distance >= min_s) && (edit_distance <= max_s)) 
                    {
                        #pragma omp critical
                        {
                            edge_lst[read_pair] = edit_distance;
                        }                    
                    }
                }
            }            
        }
    }

    // #pragma omp taskwait

}

//////////////////////
////
// void EdgeConstructor::process_blocks_in_parallel()
// {
//     int min_s = -1 * args.max_edit_dis;
//     int max_s = -1 * args.min_edit_dis;
//     // std::for_each(std::execution::par_unseq, minimiser_to_reads_.begin(), minimiser_to_reads_.end(),
//     //                 [this, min_s, max_s](const auto &entry) {
//     //                     // std::int64_t minimiser = entry.first;
//     //                     const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;
//     //                     // Process each group in parallel
//     //                     process_block(reads_vec, min_s, max_s);
//     //                 });
//     // const int num_threads = 16;  // Set the desired number of threads

//     // // Set the number of threads
//     // omp_set_num_threads(num_threads);

//     // #pragma omp parallel
//     // for (const auto &entry : minimiser_to_reads_) {
//     //     const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;
//     //     // Process each group in parallel
//     //     process_block(reads_vec, min_s, max_s);
//     // }
//     // Declare and define a global variable for available cores
//     int available_cores = omp_get_max_threads();
//     // Ensure the user-specified number of cores is within a valid range
//     int num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);

//     // Set the number of threads for OpenMP
//     omp_set_num_threads(num_cores_to_use);

//     #pragma omp parallel for
//     for (int i = 0; i < minimiser_to_reads_.size(); ++i) {
//         const auto &entry = *std::next(minimiser_to_reads_.begin(), i);
//         const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;

//         // Use OpenMP to parallelize the processing of each group of reads
//         #pragma omp task
//         {
//             // process_block(reads_vec, min_s, max_s);
//             auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s}
//                     | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_sequence1_id{}
//                     | seqan3::align_cfg::output_sequence2_id{};
//             auto alignment_results = seqan3::align_pairwise(seqan3::views::pairwise_combine(reads_vec), config);
//             auto filter_v = std::views::filter(
//                 [&min_s, &max_s](auto && res)
//                 {
//                     return (res.score() > min_s) && (res.score() <= max_s);
//                 });

//             for (auto const & result : alignment_results | filter_v)  {
//                 // seqan3::debug_stream << result << '\n';
//                 // Create a mapping, that maps from alignment-pair-ids → sequence-ids
//                 auto pair_indices = seqan3::views::pairwise_combine(std::views::iota(0ul, reads_vec.size()));
//                 // The sequence1_id is reporting the index of the alignment-pair, not the ids of the pairs them self
//                 auto alignment_id = result.sequence1_id(); // result.sequence1_id() == result.sequence2_id()

//                 auto [sid1, sid2] = pair_indices[alignment_id]; // alignment-pair-id -> sequence-id
//                 // seqan3::debug_stream << "sid1: " << sid1 << "; sid2: " << sid2 << "\n"; // the real sequence ids
//                 auto read_pair = std::make_pair(reads_vec[sid1], reads_vec[sid2]);
//                 auto it = edge_lst.find(read_pair);
//                 if (it == edge_lst.end()) {
//                     int edit_distance = result.score();
//                     #pragma omp critical
//                     {
//                         edge_lst[read_pair] = edit_distance;
//                     }                    
//                 }
//             }                
//         }            
//     }

//     #pragma omp taskwait

// }



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

// void EdgeConstructor::process_block(const std::vector<std::vector<seqan3::dna5>> &reads_vec, int min_s, int max_s)
// {
//     auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s}
//                     | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_sequence1_id{}
//                     | seqan3::align_cfg::output_sequence2_id{};

//         auto alignment_results = seqan3::align_pairwise(seqan3::views::pairwise_combine(reads_vec), config);
//         auto filter_v = std::views::filter(
//             [&min_s, &max_s](auto && res)
//             {
//                 return (res.score() > min_s) && (res.score() <= max_s);
//             });

//     // Create a mapping, that maps from alignment-pair-ids → sequence-ids
//     auto pair_indices = seqan3::views::pairwise_combine(std::views::iota(0ul, reads_vec.size()));
//     for (auto const & result : alignment_results | filter_v)  {
//         // seqan3::debug_stream << result << '\n';

//         // The sequence1_id is reporting the index of the alignment-pair, not the ids of the pairs them self
//         auto alignment_id = result.sequence1_id(); // result.sequence1_id() == result.sequence2_id()

//         auto [sid1, sid2] = pair_indices[alignment_id]; // alignment-pair-id -> sequence-id
//         // seqan3::debug_stream << "sid1: " << sid1 << "; sid2: " << sid2 << "\n"; // the real sequence ids
//         auto read_pair = std::make_pair(reads_vec[sid1], reads_vec[sid2]);
//         int edit_distance = result.score();
//         edge_lst[read_pair] = edit_distance;
//     }
//     #pragma omp critical
//     {
//         std::cout << "Processing block in thread " << omp_get_thread_num() << std::endl;
//     }
// }

// void EdgeConstructor::process_blocks_in_parallel()
// {
//     int min_s = -1 * args.max_edit_dis;
//     int max_s = -1 * args.min_edit_dis;

//     // Declare and define a global variable for available cores
//     int available_cores = omp_get_max_threads();
//     // Ensure the user-specified number of cores is within a valid range
//     int num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);

//     // Set the number of threads for OpenMP
//     omp_set_num_threads(num_cores_to_use);
//     // Identify pairs with edit distance below threshold
//     std::vector<std::pair<seqan3::dna5_vector, seqan3::dna5_vector>> remaining_read_pairs;
//     #pragma omp parallel for
//     for (auto i = 0u; i < minimiser_to_reads_.size(); ++i) {
//         const auto &entry = *std::next(minimiser_to_reads_.begin(), i);
//         const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;
//         size_t size_of_reads_vec = reads_vec.size();
//         if (size_of_reads_vec > 1000)
//         {
//             // Initialize numeric_matrix directly within the loop
//             Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> numeric_matrix;
//             size_t rows = reads_vec.size();
//             size_t cols = reads_vec[0].size();
            
//             #pragma omp parallel for
//             for (size_t i = 0; i < rows; ++i)
//             {
//                 // Convert to ranks using seqan3::views::to_rank
//                 auto rank_view = reads_vec[i] | seqan3::views::to_rank;

//                 // Transform ranks, adding 1 to each rank
//                 auto transformed_view = rank_view | views::transform([](uint8_t rank) { return rank + 1; });

//                 // Collect the result into a vector using the std::vector constructor
//                 std::vector<uint8_t> numeric_read(std::begin(transformed_view), std::end(transformed_view));

//                 // Initialize numeric_matrix if not already initialized
//                 #pragma omp critical
//                 {
//                     if (i == 0)
//                         numeric_matrix.resize(rows, numeric_read.size());
//                 }

//                 // Copy numeric_read to numeric_matrix
//                 numeric_matrix.row(i) = Eigen::Map<Eigen::Matrix<uint8_t, 1, Eigen::Dynamic>>(numeric_read.data(), 1, numeric_read.size());
//             }

//             // Calculate non-zero counts using Eigen's array-wise operations
//             // Eigen::Array<size_t, Eigen::Dynamic, Eigen::Dynamic> non_zero_counts =
//             //     (numeric_matrix.replicate(1, rows).array() != numeric_matrix.transpose().replicate(rows, 1).array()).cast<size_t>();

//             // non_zero_counts.rowwise().sum();

//             Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> result_matrix(input_matrix.rows(), input_matrix.cols());

//             // Iterate over each row
//             for (Eigen::Index i = 0; i < numeric_matrix.rows(); ++i) {
//                 // Matrix A: Replicate the current row
//                 Eigen::Matrix<int, 1, Eigen::Dynamic> matrix_A = numeric_matrix.row(i);

//                 // Matrix B: Copy the rows that come after this row
//                 Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix_B = numeric_matrix.bottomRows(numeric_matrix.rows() - (i + 1));

//                 // Calculate A - B
//                 result_matrix.row(i) = matrix_A - matrix_B.rowwise().replicate(matrix_A.cols()).transpose();

//                 // Non-zero count for the current row
//                 size_t non_zero_counts = result_matrix.array().cast<bool>().rowwise().count();

//                 #pragma omp parallel for
//                 for (Eigen::Index j = 0; j < non_zero_counts.rows(); ++j)
//                 {
//                     int count = non_zero_counts(j);
//                     if ((count > min_s) && (count <= max_s))
//                     {
//                         #pragma omp critical
//                         {
//                             edge_lst[std::make_pair(reads_vec[i], reads_vec[i + j + 1])] = count;
//                         }
//                     }
//                     else
//                     {
//                         #pragma omp critical
//                         {
//                             remaining_read_pairs.emplace_back(reads_vec[i], reads_vec[i + j + 1]);
//                         }
//                     }
//                 }

//             }

//             ////////////////////////////////////////////////////////////////////////////////
//             // std::vector<std::vector<uint8_t>> numeric_reads;
//             // #pragma omp parallel for
//             // for (auto &read : reads_vec)
//             // {
//             //     // Convert to ranks using seqan3::views::to_rank
//             //     auto rank_view = read | seqan3::views::to_rank;

//             //     // Transform ranks, adding 1 to each rank
//             //     auto transformed_view = rank_view | views::transform([](uint8_t rank) { return rank + 1; });

//             //     // Collect the result into a vector using the std::vector constructor
//             //     std::vector<uint8_t> numeric_read(std::begin(transformed_view), std::end(transformed_view));
//             //     #pragma omp critical
//             //     {
//             //         numeric_reads.push_back(numeric_read);
//             //     }
                
//             // }
//             // // #pragma omp taskwait

//             // std::vector<std::pair<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>>> remaining_read_pairs;

//             // #pragma omp parallel for
//             // for (size_t i = 0; i < numeric_reads.size(); ++i)
//             // {

//             //     // Step 2: Create matrices A and B
//             //     size_t repeat_count = numeric_reads.size() - (i + 1);
//             //     std::vector<uint8_t> repeat_element = numeric_reads[i];

//             //     // Convert repeat_element to Eigen::Matrix<uint8_t, 1, Eigen::Dynamic>
//             //     Eigen::Matrix<uint8_t, 1, Eigen::Dynamic> repeat_row =
//             //         Eigen::Map<Eigen::Matrix<uint8_t, 1, Eigen::Dynamic>>(repeat_element.data(), 1, repeat_element.size());

//             //     // Create a matrix where each row is a repeated version of repeat_row
//             //     Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> matrix_A = repeat_row.replicate(repeat_count, 1);

//             //     // Construct matrix B using index i and the remaining rows
//             //     Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> matrix_B(numeric_reads.size() - (i + 1), repeat_element.size());

//             //     for (size_t j = i + 1; j < numeric_reads.size(); ++j)
//             //     {
//             //         Eigen::Map<Eigen::Matrix<uint8_t, 1, Eigen::Dynamic>, Eigen::Unaligned> row(
//             //             numeric_reads[j].data(), 1, numeric_reads[j].size());
//             //         matrix_B.row(j - (i + 1)) = row;
//             //     }

//             //     ////////////////////////////////////
//             //     // Subtract matrices A and B
//             //     auto matrix_C = matrix_A - matrix_B;
//             //     // Step 4: Calculate non-zero counts in each row of C
//             //     std::vector<size_t> non_zero_counts(matrix_C.rows());
//             //     #pragma omp parallel for
//             //     for (Eigen::Index i = 0; i < matrix_C.rows(); ++i)
//             //     {
//             //         #pragma omp critical
//             //         {
//             //             non_zero_counts[i] = (matrix_C.row(i).array() != 0).count();
//             //         }
//             //     }

//             //     // Step 5: Determine edit distance threshold
//             //     size_t index = i + 1; // Index of the read in reads vector

//             //     // Step 6: Identify pairs with edit distance below threshold
//             //     #pragma omp parallel for
//             //     for (size_t j = 0; j < non_zero_counts.size(); ++j)
//             //     {
//             //         int count = non_zero_counts[j];
//             //         if ((count > min_s) && (count <= max_s))
//             //         {
//             //             #pragma omp critical
//             //             {
//             //                 edge_lst[std::make_pair(reads_vec[i], reads_vec[index + j])] = count;
//             //             }
//             //         }
//             //         else
//             //         {
//             //             // remaining_read_pairs.insert(std::make_pair(reads_vec[i], reads_vec[index + j]));
//             //             #pragma omp critical
//             //             {
//             //                 remaining_read_pairs.emplace_back(reads_vec[i], reads_vec[index + j]);
//             //             }
//             //         }
//             //     }
//             // }

//             // Use OpenMP to parallelize the processing of each group of reads
//             #pragma omp task
//             {
//                 // process_block(reads_vec, min_s, max_s);
//                 auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s} | seqan3::align_cfg::output_score{};

//                 #pragma omp parallel for
//                 for (auto const &pair : remaining_read_pairs)
//                 {
//                     std::vector<seqan3::dna5> seq1 = pair.first;
//                     std::vector<seqan3::dna5> seq2 = pair.second;    

//                     auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);

//                     auto read_pair = std::make_pair(seq1, seq2);
//                     auto it = edge_lst.find(read_pair);
//                     // Iterate over alignment results and access the scores
//                     for (auto const &result : alignment_results)
//                     {
//                         int edit_distance = result.score();
//                         if (it == edge_lst.end() && (edit_distance> min_s) && (edit_distance <= max_s)) 
//                         {
//                             #pragma omp critical
//                             {
//                                 edge_lst[read_pair] = edit_distance;
//                             }                    
//                         }
//                     }
//                 }            
//             }
//         }
//         else
//         {
//             // Use OpenMP to parallelize the processing of each group of reads
//             #pragma omp task
//             {
//                 // process_block(reads_vec, min_s, max_s);
//                 auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{min_s} | seqan3::align_cfg::output_score{};

//                 auto pairwise_combinations = seqan3::views::pairwise_combine(reads_vec);

//                 #pragma omp parallel for
//                 for (size_t i = 0; i < pairwise_combinations.size(); ++i)
//                 {
//                     auto const &combination = pairwise_combinations[i];
//                     auto const &seq1 = std::get<0>(combination);
//                     auto const &seq2 = std::get<1>(combination);

//                     auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);

//                     auto read_pair = std::make_pair(seq1, seq2);
//                     auto it = edge_lst.find(read_pair);
//                     // Iterate over alignment results and access the scores
//                     for (auto const &result : alignment_results)
//                     {
//                         int edit_distance = result.score();
//                         if (it == edge_lst.end() && (edit_distance> min_s) && (edit_distance <= max_s)) 
//                         {
//                             #pragma omp critical
//                             {
//                                 edge_lst[read_pair] = edit_distance;
//                             }                    
//                         }
//                     }
//                 }            
//             }
//         }
//     }
// }