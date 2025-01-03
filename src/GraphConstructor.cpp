// GraphConstructor.cpp
#include "GraphConstructor.hpp"
#include "gOMH.hpp"
#include "MinimizerGenerator.hpp"
#include "omh.hpp"
#include "StatisticsRecorder.hpp"
#include "miniception.hpp"

// #include <format>
#include <seeded_prg.hpp>
#include <boost/format.hpp>
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
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <seqan3/alphabet/all.hpp>
#include <deque>
#include <atomic>
#include <tuple>
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/depth_first_search.hpp>
using namespace boost;
// GraphConstructor::GraphConstructor(std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads, std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : key2reads_(std::move(key2reads)), read2count_(std::move(read2count)), args(args) {}
GraphConstructor::GraphConstructor(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, cmd_arguments args) : read2count_(std::move(read2count)), args(args) {}

void GraphConstructor::init_graph()
{
    // #pragma omp parallel for num_threads(args.num_process) schedule(static, 1)
    for (const auto& pair : read2count_) {
        const auto& read = pair.first;
        const auto& count = pair.second;
        auto v = boost::add_vertex({read, count}, graph_);
        read2vertex_[read] = v;
        vertex2read_[v] = read;
    }
}

void GraphConstructor::insert_edge(std::vector<seqan3::dna5> read1, std::vector<seqan3::dna5> read2, int edit_dis)
{
    auto v1 = read2vertex_[read1];
    auto v2 = read2vertex_[read2];
    if (!boost::edge(v1, v2, graph_).second) {
        // boost::add_edge(v1, v2, {read1, read2, edit_dis}, graph_);
        boost::add_edge(v1, v2, {edit_dis}, graph_);
    }
}

void GraphConstructor::edge_summary(){
    std::map<int, unsigned> weight_counts;

    // Iterate over the edges and count occurrences of each weight
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    
    for (std::tie(ei, ei_end) = boost::edges(graph_); ei != ei_end; ++ei) {
        Edge e = *ei;
        EdgeProperties& edge_props = graph_[e];
        weight_counts[edge_props.weight]++;
    }
    unsigned edge_num = 0;
    // Output the counts of edges with the same weights
    for (const auto& [weight, count] : weight_counts) {
        // std::cout << "Number of edges with weight " << weight << ": " << count << std::endl;
        Utils::getInstance().logger(LOG_LEVEL_INFO,  boost::str(boost::format("The number of edges with weight of %1%: %2%.") % weight % count));
        // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The number of edges with weight of {}: {}.", weight, count));
        edge_num += count;
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO,  boost::str(boost::format("The total number of edges: %1%.") % edge_num ));
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("The total number of edges: {}.", edge_num));
}

void GraphConstructor::construct_graph(std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads)
{
    init_graph();

    std::vector<std::vector<std::vector<seqan3::dna5>>> medium_group;
    std::vector<std::vector<std::vector<seqan3::dna5>>> large_group;
    std::atomic<int> singleton_bucket_num{0};
    auto cur_bin_n = key2reads.size();
    std::string bucket_method;
    if (args.bucketing_mode == "miniception_gomh") {
        bucket_method = "miniception";
    } else if (args.bucketing_mode == "miniception_gomh") {
        bucket_method = "random minimizer";
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets by %1%: %2%.") % bucket_method % cur_bin_n));
    // #pragma omp parallel for
    #pragma omp parallel for num_threads(args.num_process) schedule(static, 1)
    for (auto i = 0u; i < cur_bin_n; ++i) {
        const auto &entry = *std::next(key2reads.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &reads_vec = entry.second;
        auto cur_read_num = reads_vec.size();
        if ( cur_read_num == 1){
            singleton_bucket_num++;
            continue;
        } else if ( cur_read_num >= 2 && cur_read_num < args.bin_size_max){
            #pragma omp critical
            {
                // std::cout << cur_read_num << " ";
                // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("{} ", cur_read_num));
                medium_group.emplace_back(reads_vec);
            } 
        } else {
            #pragma omp critical
            {
                // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("{} ", cur_read_num));
                large_group.emplace_back(reads_vec); 
                // large_group.emplace_back(reads_vec);   
            }
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets containing one unique read after bucketing by %1%: %2%.") % bucket_method % singleton_bucket_num));
    ///////////////////////////
    // test the following function
    // auto unique_reads = mergeUniqueReads(medium_group);
    // Utils::getInstance().logger(LOG_LEVEL_INFO,  format("Test passed, unique number {}!", unique_reads.size()));

    /////////////////////////////
    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};
    // large group
    auto bucket_num = medium_group.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets with size falling in [2, %1%) by %2%: %3%.") % args.bin_size_max % bucket_method % bucket_num));
    if (bucket_num > 0){
        StatisticsRecorder recorder(bucket_num);

        int bucket_id = 0;
        for (const auto &group : medium_group)
        {
            // Atomic variables for positive and negative cases
            std::atomic<int> positive_cases{0};
            std::atomic<int> negative_cases{0};    
            auto pairwise_combinations = seqan3::views::pairwise_combine(group);
            //debug
            // std::atomic<int> total_pairs{0};
            // std::cout << "Group size: " << group.size() << std::endl;
            // std::cout << "Expected total pairs: " << (group.size() * (group.size() - 1)) / 2 << std::endl;

            // auto total_combinations = std::distance(pairwise_combinations.begin(), pairwise_combinations.end());
            // std::cout << "Total combinations calculated: " << total_combinations << std::endl;
            auto bucket_size = group.size();
            auto total_pairs = pairwise_combinations.size();
            // for (auto const &combination : pairwise_combinations) {
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (size_t i = 0; i < total_pairs; ++i){
                auto const &combination = pairwise_combinations[i];
                auto const &seq1 = std::get<0>(combination);
                auto const &seq2 = std::get<1>(combination);

                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results){
                    int edit_distance = -1 * result.score();
                    // std::cout << edit_distance << endl;
                    if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)){
                        positive_cases++;
                        #pragma omp critical
                        {
                            insert_edge(seq1, seq2, edit_distance);
                        }                    
                    } else {
                        negative_cases++;
                    }
                }
                // total_pairs++;
            }
            recorder.record_bucket(bucket_id, bucket_size, total_pairs, positive_cases, negative_cases);   
            bucket_id++;    
        } 
        // Create a timestamped filename
        std::ostringstream oss;
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm now_tm = *std::localtime(&now_time);
        if (args.bucketing_mode == "miniception_gomh") {
            oss << std::put_time(&now_tm, "miniception_bucket_stats_%Y%m%d_%H%M%S.txt");
        } else if (args.bucketing_mode == "minimizer_gomh") {
            oss << std::put_time(&now_tm, "minimizer_bucket_stats_%Y%m%d_%H%M%S.txt");
        }        
        std::filesystem::path output_file = args.output_dir / oss.str();
        recorder.write_statistics_to_file(output_file);    
        Utils::getInstance().logger(LOG_LEVEL_INFO, "Pairwise comparison for the small- or medium-size-based buckets done!");       
    } else {
        // Utils::getInstance().logger(LOG_LEVEL_INFO,  "No bucket has a size larger than 100!");
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("No bucket with size falling in [2, %1%) by %2%!") % args.bin_size_max % bucket_method));
    }
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("===== Graph construction based on %1% bucketing done! =====") % bucket_method));
    // extra large group
    if (large_group.size() > 0){
        Utils::getInstance().logger(LOG_LEVEL_INFO, "Bucketing large-size-based bins using gOMH...");
        std::mt19937_64 generator(args.seed);
        // Specify the range of values for your seeds
        std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());

        std::vector<std::vector<std::vector<seqan3::dna5>>> m_group;
        std::vector<std::vector<std::vector<seqan3::dna5>>> l_group;

        std::vector<std::pair<std::uint64_t, unsigned>> seeds_k;
        uint8_t d_t = args.max_edit_dis;
        if (args.min_edit_dis <= args.max_edit_dis){
            uint8_t times;
            if (args.max_edit_dis == 1 || args.max_edit_dis == 2){
                // times = args.omh_times - args.max_edit_dis + 1;
                times = args.gomh_times;
            } else {
                times = 1;
            }

            for (unsigned int cur_d = d_t; cur_d >= 1; cur_d--) {
                for (unsigned j = 0; j < times; ++j) {
                    std::uint64_t cur_seed = distribution(generator);
                    std::uint64_t cur_k = gOMH(args).gomh_k(args.read_length, args.probability, cur_d);
                    std::pair<std::uint64_t, unsigned> cur_pair = std::make_pair(cur_seed, cur_k);
                    seeds_k.push_back(cur_pair);                   
                } 
            }   
            /////////////////////////////////////////////////////
            // When the permutation_times larger than the number of k-mer candidates and the kmer size are the same one, bucketing the reads using each kmer candidate
            int flag = 0;
            unsigned first_k = seeds_k.front().second;
            for (const auto& pair : seeds_k) {
                if (pair.second != first_k) {
                    flag += 1;
                }
            }
            if (flag == 0) {
                auto gomh_kmer_n = static_cast<int>(args.read_length - 2 * first_k + 1);
                auto permutation_times = args.max_edit_dis * times;  
                if (permutation_times >= gomh_kmer_n){
                    args.gomh_flag = true;
                }              
            }
            /////////////////////////////////////////////////////
        } else {
             Utils::getInstance().logger(LOG_LEVEL_ERROR, boost::str(boost::format("min_edit_dis(%1%) should not be larger than max_edit_dis(%2%)") % args.min_edit_dis % args.max_edit_dis));
        }
        std::atomic<int> singleton_bucket_num{0};
        for (const auto &el_group : large_group){
            auto cur_hash2reads = gOMH(args).gomh2read_main(el_group, seeds_k);
            auto cur_bin_n = cur_hash2reads.size();
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (auto i = 0u; i < cur_bin_n; ++i) {
                const auto &cur_entry = *std::next(cur_hash2reads.begin(), i);
                const std::vector<std::vector<seqan3::dna5>> &cur_reads_vec = cur_entry.second;
                auto cur_num = cur_reads_vec.size();
                if ( cur_num == 1){
                    singleton_bucket_num++;
                    continue;
                } else if (cur_bin_n >= 2 && cur_bin_n < args.bin_size_max){
                    #pragma omp critical
                    {
                        // std::cout << cur_num << " ";
                        // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("{} ", cur_num));
                        m_group.emplace_back(cur_reads_vec);
                    } 
                } else {
                    #pragma omp critical
                    l_group.emplace_back(cur_reads_vec);
                } 
            }
        }
        args.gomh_flag = false; // make this flag false after using gomh2read_main
        auto bucket_num = m_group.size();
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("After bucketing large-size groups with gOMH, there are %1% single-size buckets and %2% buckets ranging from size 2 to %3%.") % singleton_bucket_num % bucket_num % args.bin_size_max));
        if (bucket_num > 0){
            StatisticsRecorder recorder(bucket_num);
            int bucket_id = 0;
            for (const auto &cur_reads_vec : m_group){
                // Atomic variables for positive and negative cases
                std::atomic<int> positive_cases{0};
                std::atomic<int> negative_cases{0};  
                auto pairwise_combinations = seqan3::views::pairwise_combine(cur_reads_vec);
                auto total_pairs = pairwise_combinations.size();
                // #pragma omp parallel for
                #pragma omp parallel for num_threads(args.num_process) schedule(static)
                for (size_t i = 0; i < total_pairs; ++i)
                {
                    auto const &combination = pairwise_combinations[i];
                    auto const &seq1 = std::get<0>(combination);
                    auto const &seq2 = std::get<1>(combination);
                    auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                    // Iterate over alignment results and access the scores
                    for (auto const &result : alignment_results)
                    {
                        int edit_distance = -1 * result.score();
                        // std::cout << edit_distance << endl;
                        //if ((edit_distance >= min_s) && (edit_distance <= max_s))
                        if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)) 
                        // if ((edit_distance < min_s) && (edit_distance > max_s))
                        {
                            positive_cases++;
                            #pragma omp critical
                            {
                                // edge_lst[read_pair_set] = edit_distance;
                                insert_edge(seq1, seq2, edit_distance);
                            }                    
                        } else {
                            negative_cases++;
                        }
                    }
                }
                auto bucket_size = cur_reads_vec.size();
                recorder.record_bucket(bucket_id, bucket_size, total_pairs, positive_cases, negative_cases);   
                bucket_id++;       
            }
            std::ostringstream oss;
            auto now = std::chrono::system_clock::now();
            std::time_t now_time = std::chrono::system_clock::to_time_t(now);
            std::tm now_tm = *std::localtime(&now_time);
            oss << std::put_time(&now_tm, "gomh_bucket_stats_%Y%m%d_%H%M%S.txt");
            std::filesystem::path output_file = args.output_dir / oss.str();
            recorder.write_statistics_to_file(output_file);             
            // Utils::getInstance().logger(LOG_LEVEL_INFO,  "Graph update for the small- or medium-size-based buckets generated by gOMH bucketing done!"); 
            edge_summary(); 
            Utils::getInstance().logger(LOG_LEVEL_INFO, "===== Graph update based on gOMH bucketing done ====="); 
        }
    
        //////////////////////////////////////////////////////////////////////////////////////////////
        if (l_group.size() > 0){
            auto unique_reads = mergeUniqueReads(l_group);
            Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of remaining unprocessed unique reads: %1%") % unique_reads.size()));
            //////////////////////////
            // Method 1
            std::vector<std::pair<int, int>> v_pairs;
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (const auto &cur_read : unique_reads){
                auto cur_v = read2vertex_[cur_read];
                auto indirect_neighbors = visitNeighborsOfNeighborsWithThreshold(graph_, cur_v, args.visit_depth);
                if (!indirect_neighbors.empty()) {
                    for (auto v : indirect_neighbors){
                        std::pair<int, int> cur_pair = std::make_pair(cur_v, v);
                        #pragma omp critical
                        v_pairs.emplace_back(cur_pair);
                    }
                    // std::cout << "good" << endl;
                }
            }
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (const auto& pair : v_pairs) {
                auto seq1 = vertex2read_[pair.first];
                auto seq2 = vertex2read_[pair.second];
                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results)
                {
                    int edit_distance = -1 * result.score();
                    //if ((edit_distance >= min_s) && (edit_distance <= max_s))
                    if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)) 
                    {
                        #pragma omp critical
                        {
                            insert_edge(seq1, seq2, edit_distance);
                        }                    
                    } 
                } 
            }
            edge_summary(); 
            Utils::getInstance().logger(LOG_LEVEL_INFO, "===== Graph update for the remaining unprocessed unique reads with graph traversal done! =====");
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////
    if ((args.read_length < 16) && (medium_group.size() > 0)){
        auto medium_unique_reads = mergeUniqueReads(medium_group);
        update_graph_omh(medium_unique_reads); 
        Utils::getInstance().logger(LOG_LEVEL_INFO, "===== Graph update with additional gOMH done! =====");
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // if (args.min_edit_dis > 1 && args.read_length >= 50) {
    if (args.min_edit_dis > 1) {
        remove_edges_in_interval(graph_, 1, args.min_edit_dis - 1);
        edge_summary();
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("===== Graph update for removing the edges with weights less than %1%. =====") % args.min_edit_dis));
    }
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Edit-distance-based read graph construction done!");
    Utils::getInstance().logger(LOG_LEVEL_INFO, "*************************************************");
}

void GraphConstructor::update_graph_omh(std::vector<std::vector<seqan3::dna5>> unique_reads){
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Additional bucketing of reads in medium-sized buckets produced by minimizer or miniception using gOMH...!");
    std::mt19937_64 generator(args.seed + 1);
    // Specify the range of values for your seeds
    std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());

    std::vector<std::vector<std::vector<seqan3::dna5>>> m_group;
    std::vector<std::vector<std::vector<seqan3::dna5>>> l_group;

    std::vector<std::pair<std::uint64_t, unsigned>> seeds_k;
    uint8_t d_t = args.max_edit_dis;
    if (args.min_edit_dis <= args.max_edit_dis){
        uint8_t times;
        if (args.max_edit_dis == 1 || args.max_edit_dis == 2){
            // times = args.omh_times - args.max_edit_dis + 1;
            times = args.gomh_times;
        } else {
            times = 1;
        }
        for (unsigned int cur_d = d_t; cur_d >= 1; cur_d--) {
            for (unsigned j = 0; j < times; ++j) {
                std::uint64_t cur_seed = distribution(generator);
                std::uint64_t cur_k = gOMH(args).gomh_k(args.read_length, args.probability, cur_d);
                std::pair<std::uint64_t, unsigned> cur_pair = std::make_pair(cur_seed, cur_k);
                seeds_k.push_back(cur_pair);                   
            } 
        }   

        /////////////////////////////////////////////////////
        // When the permutation_times larger than the number of k-mer candidates and the kmer size are the same one, bucketing the reads using each kmer candidate
        int flag = 0;
        unsigned first_k = seeds_k.front().second;
        for (const auto& pair : seeds_k) {
            if (pair.second != first_k) {
                flag += 1;
            }
        }
        if (flag == 0) {
            auto gomh_kmer_n = static_cast<int>(args.read_length - 2 * first_k + 1);
            auto permutation_times = args.max_edit_dis * times;  
            if (permutation_times >= gomh_kmer_n){
                args.gomh_flag = true;
            }              
        }
        /////////////////////////////////////////////////////

    } else {
        // Utils::getInstance().logger(LOG_LEVEL_ERROR,  std::format("min_edit_dis({}) should not larger than max_edit_dis({}) ", args.min_edit_dis, args.max_edit_dis));
        Utils::getInstance().logger(LOG_LEVEL_ERROR, boost::str(boost::format("min_edit_dis(%1%) should not be larger than max_edit_dis(%2%)") % args.min_edit_dis % args.max_edit_dis));

    }

    auto cur_hash2reads = gOMH(args).gomh2read_main(unique_reads, seeds_k);
    args.gomh_flag = false; // make this flag false after using gomh2read_main
    auto cur_bin_n = cur_hash2reads.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets generated by gOMH is %1%.") % cur_bin_n));
    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};
    std::atomic<int> singleton_bucket_num{0};
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto i = 0u; i < cur_bin_n; ++i) {
        const auto &cur_entry = *std::next(cur_hash2reads.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &cur_reads_vec = cur_entry.second;
        auto cur_num = cur_reads_vec.size();
        if ( cur_num == 1){
            singleton_bucket_num++;
            continue;
        } else if (cur_bin_n >= 2 && cur_bin_n < args.bin_size_max){
            #pragma omp critical
            {
                m_group.emplace_back(cur_reads_vec);
            } 
        } else {
            #pragma omp critical
            l_group.emplace_back(cur_reads_vec);
        } 
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets containing one unique read after bucketing by gOMH is %1%.") % singleton_bucket_num));
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Bucketing unique reads using gOMH done!");
    auto bucket_num = m_group.size();
    if (bucket_num > 0){
        StatisticsRecorder recorder(bucket_num);
        int bucket_id = 0;
        for (const auto &cur_reads_vec : m_group){
            // Atomic variables for positive and negative cases
            std::atomic<int> positive_cases{0};
            std::atomic<int> negative_cases{0};   
            auto pairwise_combinations = seqan3::views::pairwise_combine(cur_reads_vec);
            auto bucket_size = cur_reads_vec.size();
            auto total_pairs = pairwise_combinations.size();
            // #pragma omp parallel for
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (size_t i = 0; i < total_pairs; ++i)
            {
                auto const &combination = pairwise_combinations[i];
                auto const &seq1 = std::get<0>(combination);
                auto const &seq2 = std::get<1>(combination);
                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results)
                {
                    int edit_distance = -1 * result.score();
                    if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)){
                        positive_cases++;
                        #pragma omp critical
                        {
                            insert_edge(seq1, seq2, edit_distance);
                        }                    
                    } else {
                        negative_cases++;
                    }
                }
            }
            recorder.record_bucket(bucket_id, bucket_size, total_pairs, positive_cases, negative_cases);   
            bucket_id++;     
        }
        // Create a timestamped filename
        std::ostringstream oss;
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::tm now_tm = *std::localtime(&now_time);
        oss << std::put_time(&now_tm, "additional_gomh_bucket_stats_%Y%m%d_%H%M%S.txt");
        std::filesystem::path output_file = args.output_dir / oss.str();
        recorder.write_statistics_to_file(output_file);      
        edge_summary();   
        Utils::getInstance().logger(LOG_LEVEL_INFO, "Graph update for the small- or medium-size-based buckets generated by gOMH bucketing done!");   
    }
    // 

    //////////////////////////////////////////////////////////////////////////////////////////////
    if (l_group.size() > 0){
        auto unique_reads = mergeUniqueReads(l_group);
        // Utils::getInstance().logger(LOG_LEVEL_INFO,   std::format("The number of remaining unprocessed unique reads: {} ", unique_reads.size()));
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of remaining unprocessed unique reads: %1%") % unique_reads.size()));
        //////////////////////////
        // Method 1
        std::vector<std::pair<int, int>> v_pairs;
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (const auto &cur_read : unique_reads){
            auto cur_v = read2vertex_[cur_read];
            auto indirect_neighbors = visitNeighborsOfNeighborsWithThreshold(graph_, cur_v, args.visit_depth);
            if (!indirect_neighbors.empty()) {
                for (auto v : indirect_neighbors){
                    std::pair<int, int> cur_pair = std::make_pair(cur_v, v);
                    #pragma omp critical
                    v_pairs.emplace_back(cur_pair);
                }
                // std::cout << "good" << endl;
            }
        }
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (const auto& pair : v_pairs) {
            auto seq1 = vertex2read_[pair.first];
            auto seq2 = vertex2read_[pair.second];
            auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
            // Iterate over alignment results and access the scores
            for (auto const &result : alignment_results)
            {
                int edit_distance = -1 * result.score();
                if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)) 
                {
                    #pragma omp critical
                    {
                        insert_edge(seq1, seq2, edit_distance);
                    }                    
                } 
            } 
        } 
        edge_summary();
        Utils::getInstance().logger(LOG_LEVEL_INFO, "Graph update for the remaining unprocessed unique reads done!");
    }
}

std::vector<std::vector<seqan3::dna5>> GraphConstructor::mergeUniqueReads(const std::vector<std::vector<std::vector<seqan3::dna5>>>& read_vectors) {
    if (read_vectors.size() == 1) {
        return read_vectors[0];
    } else {
        std::vector<std::vector<seqan3::dna5>> merged_reads;
        size_t total_reads_count = 0;
        for (const auto& read_vector : read_vectors) {
            total_reads_count += read_vector.size();
        }
        merged_reads.reserve(total_reads_count);

        // Merge all read vectors into the merged_reads vector in parallel
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (size_t i = 0; i < read_vectors.size(); ++i) {
                    #pragma omp task
                    {
                        const auto& read_vector = read_vectors[i];
                        for (const auto& read : read_vector) {
                            #pragma omp critical
                            merged_reads.push_back(read);
                        }
                    }
                }
            }
        }

        // Sort the merged vector of reads in parallel
        #pragma omp parallel
        {
            #pragma omp single nowait
            std::sort(merged_reads.begin(), merged_reads.end());
        }

        // Remove duplicates using std::unique and erase idiom in parallel
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                auto last = std::unique(merged_reads.begin(), merged_reads.end());
                #pragma omp task
                merged_reads.erase(last, merged_reads.end());
            }
        }

        return merged_reads;
    }
}

// Function to visit neighbors of a given node until distance exceeds threshold
void GraphConstructor::visitNeighborsWithThreshold(const Graph& g, Vertex node, int distance_threshold, int current_distance, std::vector<Vertex>& indirect_neighbors, std::vector<bool>& visited) {
    // Mark the current node as visited
    visited[node] = true;

    // Iterate over the adjacent vertices of the given node
    graph_traits<Graph>::adjacency_iterator ai, ai_end;
    for (boost::tie(ai, ai_end) = adjacent_vertices(node, g); ai != ai_end; ++ai) {
        Vertex neighbor = *ai;
        // Visit neighbor if not visited and distance does not exceed threshold
        if (!visited[neighbor] && current_distance + 1 <= distance_threshold) {
            indirect_neighbors.push_back(neighbor);
            // Recursively visit neighbors of neighbors
            visitNeighborsWithThreshold(g, neighbor, distance_threshold, current_distance + 1, 
                                         indirect_neighbors, visited);
        } else {
            return;
        }
    }
}

// Function to visit neighbors of neighbors until distance exceeds a threshold
std::vector<Vertex> GraphConstructor::visitNeighborsOfNeighborsWithThreshold(const Graph& g, Vertex node, int distance_threshold) {
    std::vector<Vertex> indirect_neighbors;
    std::vector<bool> visited(num_vertices(g), false); // Initialize visited array

    // Visit neighbors of neighbors with the specified threshold
    visitNeighborsWithThreshold(g, node, distance_threshold, 0, indirect_neighbors, visited);

    return indirect_neighbors;
}

void GraphConstructor::save_graph() const {
    // std::ofstream out(graph_full_path_);
    // auto graph_full_path_ = args.output_dir / args.graph_filename;
    std::filesystem::path modifiable_path = args.input_data;
    modifiable_path.replace_extension(".dot");
    std::string output_filename = modifiable_path.filename().string();
    auto graph_full_path_ = args.output_dir / output_filename;
    std::cout << "Output graph Path: " << graph_full_path_ << std::endl;
    // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Output graph Path: {}.", graph_full_path_));
    std::ofstream dot_file(graph_full_path_);

    // boost::write_graphviz(dot_file, graph);  
    
    // Check if the file is open
    if (dot_file.is_open()) {
        // Write the graph to the file in Graphviz format
        // boost::write_graphviz(dot_file, graph, make_label_writer(get(&VertexProperties::count, graph)),
        //                make_label_writer(get(&EdgeProperties::weight, graph)));

        // boost::write_graphviz(dot_file, graph, make_label_writer(get(&VertexProperties::count, graph)), custom_edge_label_writer<Graph>(graph));
        boost::write_graphviz(dot_file, graph_, custom_vertex_label_writer<Graph>(graph_), custom_edge_label_writer<Graph>(graph_));

        // Close the file
        dot_file.close();
        
        std::cout << "Graph written to " << graph_full_path_ << " successfully." << std::endl;
        // Utils::getInstance().logger(LOG_LEVEL_DEBUG,  std::format("Graph written to {} successfully.", graph_full_path_));
    } else {
        std::cerr << "Error opening file " << graph_full_path_ << std::endl;
    }        
    // std::cout << "Graph saved to file successfully." << std::endl;
}

void GraphConstructor::construt_graph_via_pairwise_comparison(std::vector<std::vector<seqan3::dna5>> unique_reads)
{
    init_graph();

    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};

    auto pairwise_combinations = seqan3::views::pairwise_combine(unique_reads);

    #pragma omp parallel for num_threads(args.num_process) schedule(static, 1)
    for (size_t i = 0; i < pairwise_combinations.size(); ++i)
    {
        auto const &combination = pairwise_combinations[i];
        auto const &seq1 = std::get<0>(combination);
        auto const &seq2 = std::get<1>(combination);

        auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
        // Iterate over alignment results and access the scores
        for (auto const &result : alignment_results)
        {
            int edit_distance = -1 * result.score();
            if ((edit_distance >= args.min_edit_dis) && (edit_distance <= args.max_edit_dis))
            {
                #pragma omp critical
                {
                    insert_edge(seq1, seq2, edit_distance); 
                }                    
            }
        }
    }
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Constructing edit-distance-based read graph via brute force done!");
    Utils::getInstance().logger(LOG_LEVEL_INFO, "*****************************************************************");
}

void GraphConstructor::construt_graph_via_miniception(std::vector<std::vector<seqan3::dna5>> unique_reads)
{
    init_graph();

    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};
    
    auto cur_hash2reads = Miniception(args).miniception2read_main(unique_reads, args.miniception_k, args.miniception_w,  args.num_process);

    auto cur_bin_n = cur_hash2reads.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets: %1%.") % cur_bin_n));

    StatisticsRecorder recorder(cur_bin_n);
    int bucket_id = 0;
    int singleton_bucket_num = 0;
    for (auto i = 0u; i < cur_bin_n; ++i) {
        const auto &cur_entry = *std::next(cur_hash2reads.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &cur_reads_vec = cur_entry.second;
        auto cur_num = cur_reads_vec.size();
        if ( cur_num == 1){
            singleton_bucket_num++;
            continue;
        } else if (cur_num >= 2){ 
            // Atomic variables for positive and negative cases
            std::atomic<int> positive_cases{0};
            std::atomic<int> negative_cases{0};  
            auto pairwise_combinations = seqan3::views::pairwise_combine(cur_reads_vec);
            auto total_pairs = pairwise_combinations.size();
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (size_t i = 0; i < total_pairs; ++i)
            {
                auto const &combination = pairwise_combinations[i];
                auto const &seq1 = std::get<0>(combination);
                auto const &seq2 = std::get<1>(combination);
                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results)
                {
                    int edit_distance = -1 * result.score();
                    if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)) 
                    {
                        positive_cases++;
                        #pragma omp critical
                        {
                            // edge_lst[read_pair_set] = edit_distance;
                            insert_edge(seq1, seq2, edit_distance);
                        }                    
                    } else {
                        negative_cases++;
                    }
                }
            }  
            recorder.record_bucket(bucket_id, cur_num, total_pairs, positive_cases, negative_cases);   
            bucket_id++;    
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets containing one unique read: %1%.") % singleton_bucket_num));
    std::ostringstream oss;
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_time);
    oss << std::put_time(&now_tm, "miniception_bucket_stats_%Y%m%d_%H%M%S.txt");
    std::filesystem::path output_file = args.output_dir / oss.str();
    recorder.write_statistics_to_file(output_file);  
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Constructing edit-distance-based read graph via miniception bucketing done!");
    Utils::getInstance().logger(LOG_LEVEL_INFO, "******************************************************************************");
}

void GraphConstructor::construt_graph_via_omh(std::vector<std::vector<seqan3::dna5>> unique_reads)
{
    init_graph();

    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};

    std::mt19937_64 generator(args.seed);
    std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());
    std::uint64_t cur_seed = distribution(generator);
    omh_sketcher<std::mt19937_64> sketcher(args.omh_k, args.omh_l, cur_seed);

    // Group reads by OMH sketches
    OMH omh;
    std::mt19937_64 prg(args.seed);
    auto cur_hash2reads = omh.omh2read_main(unique_reads, sketcher, prg, args.omh_m, args.num_process);
    auto cur_bin_n = cur_hash2reads.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets: %1%.") % cur_bin_n));

    StatisticsRecorder recorder(cur_bin_n);
    int bucket_id = 0;
    int singleton_bucket_num = 0;
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto i = 0u; i < cur_bin_n; ++i) {
        const auto &cur_entry = *std::next(cur_hash2reads.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &cur_reads_vec = cur_entry.second;
        auto cur_num = cur_reads_vec.size();
        if ( cur_num == 1){
            singleton_bucket_num++;
            continue;
        } else if (cur_num >= 2){ 
            // Atomic variables for positive and negative cases
            std::atomic<int> positive_cases{0};
            std::atomic<int> negative_cases{0};  

            auto pairwise_combinations = seqan3::views::pairwise_combine(cur_reads_vec);
            auto total_pairs = pairwise_combinations.size();
            // #pragma omp parallel for
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (size_t i = 0; i < total_pairs; ++i)
            {
                auto const &combination = pairwise_combinations[i];
                auto const &seq1 = std::get<0>(combination);
                auto const &seq2 = std::get<1>(combination);
                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results)
                {
                    int edit_distance = -1 * result.score();
                    if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)) 
                    {
                        positive_cases++;
                        #pragma omp critical
                        {
                            insert_edge(seq1, seq2, edit_distance);
                        }                    
                    } else {
                        negative_cases++;
                    }
                }
            }  
            recorder.record_bucket(bucket_id, cur_num, total_pairs, positive_cases, negative_cases);   
            bucket_id++;    
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets containing one unique read is %1%.") % singleton_bucket_num));

    std::ostringstream oss;
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_time);
    oss << std::put_time(&now_tm, "omh_bucket_stats_%Y%m%d_%H%M%S.txt");
    std::filesystem::path output_file = args.output_dir / oss.str();
    recorder.write_statistics_to_file(output_file);   
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Constructing edit-distance-based read graph via original OMH done!");
    Utils::getInstance().logger(LOG_LEVEL_INFO, "******************************************************************");
}

/* 
void GraphConstructor::construt_graph_via_minimizer_only(std::vector<std::vector<seqan3::dna5>> unique_reads)
{
    init_graph();

    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};

    MinimizerGenerator minimizer_generator(args);
    auto cur_hash2reads = minimizer_generator.minimizer_only2reads_main(unique_reads, args.minimizer_k, args.minimizer_m);

    auto cur_bin_n = cur_hash2reads.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets: %1%.") % cur_bin_n));

    StatisticsRecorder recorder(cur_bin_n);
    int bucket_id = 0;
    int singleton_bucket_num = 0;
    for (auto i = 0u; i < cur_bin_n; ++i) {
        const auto &cur_entry = *std::next(cur_hash2reads.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &cur_reads_vec = cur_entry.second;
        auto cur_num = cur_reads_vec.size();
        if ( cur_num == 1){
            singleton_bucket_num++;
            continue;
        } else if (cur_num >= 2){ 
            // Atomic variables for positive and negative cases
            std::atomic<int> positive_cases{0};
            std::atomic<int> negative_cases{0};  
            auto pairwise_combinations = seqan3::views::pairwise_combine(cur_reads_vec);
            auto total_pairs = pairwise_combinations.size();
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for (size_t i = 0; i < total_pairs; ++i)
            {
                auto const &combination = pairwise_combinations[i];
                auto const &seq1 = std::get<0>(combination);
                auto const &seq2 = std::get<1>(combination);
                auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);
                // Iterate over alignment results and access the scores
                for (auto const &result : alignment_results)
                {
                    int edit_distance = -1 * result.score();
                    if ((edit_distance >= 1) && (edit_distance <= args.max_edit_dis)) 
                    {
                        positive_cases++;
                        #pragma omp critical
                        {
                            // edge_lst[read_pair_set] = edit_distance;
                            insert_edge(seq1, seq2, edit_distance);
                        }                    
                    } else {
                        negative_cases++;
                    }
                }
            }  
            recorder.record_bucket(bucket_id, cur_num, total_pairs, positive_cases, negative_cases);   
            bucket_id++;    
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets containing one unique read: %1%.") % singleton_bucket_num));
    std::ostringstream oss;
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_time);
    oss << std::put_time(&now_tm, "minimizer_only_bucket_stats_%Y%m%d_%H%M%S.txt");
    std::filesystem::path output_file = args.output_dir / oss.str();
    recorder.write_statistics_to_file(output_file);  
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Constructing edit-distance-based read graph via minimizer bucketing only done!");
    Utils::getInstance().logger(LOG_LEVEL_INFO, "******************************************************************************");
}
*/

/*
void GraphConstructor::construt_graph_via_omh(std::vector<std::vector<seqan3::dna5>> unique_reads)
{
    init_graph();

    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-1 * args.max_edit_dis} | seqan3::align_cfg::output_score{};

    std::mt19937_64 generator(args.seed);
    std::uniform_int_distribution<std::uint64_t> distribution(std::numeric_limits<std::uint64_t>::min(), std::numeric_limits<std::uint64_t>::max());
    std::uint64_t cur_seed = distribution(generator);
    omh_sketcher<std::mt19937_64> sketcher(args.omh_k, args.omh_l, cur_seed);

    // Group reads by OMH sketches
    OMH omh;
    std::mt19937_64 prg(args.seed);
    auto cur_hash2reads = omh.omh2read_main(unique_reads, sketcher, prg, args.omh_m, args.num_process);
    auto cur_bin_n = cur_hash2reads.size();
    Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets: %1%.") % cur_bin_n));
    // StatisticsRecorder recorder(cur_bin_n);
    // int bucket_id = 0;
    // int singleton_bucket_num = 0;
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto i = 0u; i < cur_bin_n; ++i) {
        const auto &cur_entry = *std::next(cur_hash2reads.begin(), i);
        const std::vector<std::vector<seqan3::dna5>> &cur_reads_vec = cur_entry.second;
        auto cur_num = cur_reads_vec.size();
        if ( cur_num == 1){
            continue;
        } else if (cur_num >= 2){ 
            // Atomic variables for positive and negative cases
            std::atomic<int> positive_cases{0};
            std::atomic<int> negative_cases{0};  

            auto pairwise_combinations = seqan3::views::pairwise_combine(cur_reads_vec);
            auto total_pairs = pairwise_combinations.size();

            #pragma omp parallel num_threads(args.num_process) schedule(dynamic)
            {
                // Thread-local storage for edges to insert
                std::vector<std::tuple<std::vector<seqan3::dna5>, std::vector<seqan3::dna5>, int>> local_edges;

                #pragma omp for schedule(dynamic)
                for (size_t i = 0; i < total_pairs; ++i) {
                    auto const& combination = pairwise_combinations[i];
                    auto const& seq1 = std::get<0>(combination);
                    auto const& seq2 = std::get<1>(combination);
                    auto alignment_results = seqan3::align_pairwise(std::tie(seq1, seq2), config);

                    for (auto const& result : alignment_results) {
                        int edit_distance = -1 * result.score();
                        if (edit_distance >= 1 && edit_distance <= args.max_edit_dis) {
                            local_edges.emplace_back(seq1, seq2, edit_distance);
                        }
                    }
                }

                // Critical section for merging local edges into global structure
                #pragma omp critical
                {
                    for (auto const& edge : local_edges) {
                        const auto& [read1, read2, edit_dis] = edge;
                        insert_edge(read1, read2, edit_dis);
                    }
                }
            }  
        }
    }
    // Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of buckets containing one unique read is %1%.") % singleton_bucket_num));

    // std::ostringstream oss;
    // auto now = std::chrono::system_clock::now();
    // std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    // std::tm now_tm = *std::localtime(&now_time);
    // oss << std::put_time(&now_tm, "omh_bucket_stats_%Y%m%d_%H%M%S.txt");
    // std::filesystem::path output_file = args.output_dir / oss.str();
    // recorder.write_statistics_to_file(output_file);   
    edge_summary();
    Utils::getInstance().logger(LOG_LEVEL_INFO, "Constructing edit-distance-based read graph via original OMH done!");
    Utils::getInstance().logger(LOG_LEVEL_INFO, "******************************************************************");
}

*/