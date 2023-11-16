/**
 * @brief Construct a new Read Write:: Read Write object
 * 
 */
#include "ReadWrite.h"

ReadWrite::ReadWrite(void){

}

ReadWrite::~ReadWrite(void){

}

std::pair<std::set<std::vector<seqan3::dna5>>, std::map<std::vector<seqan3::dna5>, uint32_t>> ReadWrite::get_unique_reads_counts(cmd_arguments args){
    seqan3::sequence_file_input fin{args.input_data};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};

    // Define a set to store unique reads.
    std::set<std::vector<seqan3::dna5>> unique_reads;
    // Define a map to store unique reads and their counts
    std::map<std::vector<seqan3::dna5>, uint32_t> unique_reads2counts;

    // You can use a for loop:
    for (auto & record : fin)
    {
        // auto seq = record.sequence() | seqan3::views::to_char;
        // // seqan3::debug_stream << seq << '\n';
        // string seq_str(seq.begin(), seq.end());
        // // std::cout << seq_str << endl;
        // unique_reads_counts[seq_str]++;
        unique_reads2counts[record.sequence()]++;
        // Insert the sequence into the set to keep it unique.
        unique_reads.insert(record.sequence());
    }
    // return unique_reads; 
    return {unique_reads, unique_reads2counts};   
}


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
