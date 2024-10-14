std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH::ori_omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k){
    if (args.gomh_flag){
        // When the permutation_times larger than the number of k-mer candidates and the kmer size are the same one, bucketing the reads using each kmer candidate
        auto first_pair = seeds_k[0];
        std::uint64_t seed = first_pair.first;
        unsigned k = first_pair.second;
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (auto const & read : unique_reads){  
            auto gomh_values = omh_pos2(read, k, seed);
            for (auto const & gomh_val : gomh_values){
                #pragma omp critical
                {
                    gomh2reads[gomh_val].push_back(read);
                }                    
            }
        }  
    } else {
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (auto const & read : unique_reads){
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for(auto &pair : seeds_k){
                std::uint64_t seed = pair.first;
                unsigned k = pair.second;
                auto omh_value = ori_omh_pos(read, k, seed); 
                #pragma omp critical
                {
                    omh2reads[omh_value].push_back(read);
                }        
            } 
        }        
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
    return omh2reads;            
}