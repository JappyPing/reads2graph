/**
 * @brief Construct a new Read Write:: Read Write object
 * 
 */
#include "ReadWrite.hpp"
#include "LoggingLevels.hpp"

// #include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <seqan3/io/sequence_file/all.hpp>

ReadWrite::ReadWrite(cmd_arguments args) : args(args){
}

ReadWrite::~ReadWrite(void){

}

// std::pair<std::vector<std::vector<seqan3::dna5>>, std::unordered_map<std::vector<seqan3::dna5>, uint32_t>> ReadWrite::get_unique_reads_counts(){
//     Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Input dataset: {} ", args.input_data.string()));

//     seqan3::sequence_file_input file_in{args.input_data};
//     std::vector<decltype(file_in)::record_type> records;
//     for (auto &record : file_in)
//     {
//         records.push_back(std::move(record));
//     }

//     std::unordered_map<std::vector<seqan3::dna5>, uint32_t> read2count;
//     std::vector<std::vector<seqan3::dna5>> unique_reads;

//     #pragma omp parallel for
//     for (auto &record : records)
//     {
//         // Process the record and update local data structure
//         std::vector<seqan3::dna5> cur_seq = record.sequence();

//         // Combine local result into global data structure using critical section
//         #pragma omp critical
//         {
//             if (read2count[cur_seq]++ == 0) {
//                 unique_reads.push_back(cur_seq);
//             }
//         }
//     }
//     return {unique_reads, read2count};   
// }
std::pair<std::vector<std::vector<seqan3::dna5>>, std::map<std::vector<seqan3::dna5>, uint32_t>> ReadWrite::get_unique_reads_counts(){
// std::map<std::vector<seqan3::dna5>, uint32_t> ReadWrite::get_unique_reads_counts(){
    Utils::getInstance().logger(LOG_LEVEL_INFO,  std::format("Input dataset: {} ", args.input_data.string()));
    // std::cout << "Input dataset: " << args.input_data.string() << endl;
    seqan3::sequence_file_input fin{args.input_data};
    // using record_type = decltype(fin)::record_type;
    // std::vector<record_type> records{};
    std::vector<decltype(fin)::record_type> records;
    for (auto &record : fin)
    {
        records.push_back(std::move(record));
    }
    // Define a set to store unique reads.
    std::vector<std::vector<seqan3::dna5>> unique_reads;
    // Define a map to store unique reads and their counts
    std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    // std::unordered_map<seqan3::dna5_vector, std::string> read2id;
    
    // OpenMP parallel for loop
    // OpenMP parallel for loop using range-based for loop
    #pragma omp parallel for
    for (auto &record : records)
    {
        // Process the record and update local data structure
        std::vector<seqan3::dna5> cur_seq = record.sequence();

        // Combine local result into global data structure using critical section
        #pragma omp critical
        {
            // Combine unique_reads
            // unique_reads.insert(cur_seq);
            if (read2count[cur_seq]++ == 0) {
                unique_reads.push_back(cur_seq);
            }
            // Combine read2count
            // read2count[current_sequence]++;
        }
    }
    // return unique_reads; 
    return {unique_reads, read2count};   
    // return read2count; 
}


// std::pair<std::set<std::vector<seqan3::dna5>>, std::map<std::vector<seqan3::dna5>, uint32_t>> ReadWrite::get_unique_reads_counts(cmd_arguments args){
//     seqan3::sequence_file_input fin{args.input_data};
//     // using record_type = decltype(fin)::record_type;
//     // std::vector<record_type> records{};

//     // Define a set to store unique reads.
//     std::set<std::vector<seqan3::dna5>> unique_reads;
//     // Define a map to store unique reads and their counts
//     std::map<std::vector<seqan3::dna5>, uint32_t> unique_reads2counts;

//     // You can use a for loop:
//     for (auto & record : fin)
//     {
//         // auto seq = record.sequence() | seqan3::views::to_char;
//         // // seqan3::debug_stream << seq << '\n';
//         // string seq_str(seq.begin(), seq.end());
//         // // std::cout << seq_str << endl;
//         // unique_reads_counts[seq_str]++;
//         unique_reads2counts[record.sequence()]++;
//         // Insert the sequence into the set to keep it unique.
//         unique_reads.insert(record.sequence());
//     }
//     // return unique_reads; 
//     return {unique_reads, unique_reads2counts};   
// }


    // // // Replace with the path to your output file for unique reads and counts
    // std::filesystem::path output_file_path = "../data/unique_reads_counts.txt";

    // // Open the output file for writing
    // std::ofstream output_file(output_file_path);
    // // Write unique reads and their counts to the output file
    // for (auto const & [sequence, count] : unique_reads_counts)
    // {
    //     output_file << sequence << '\t' << count << '\n';
    // }

    // seqan3::debug_stream << "Unique reads and their counts written to " << output_file_path << '\n';


