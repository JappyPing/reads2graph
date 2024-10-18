#ifndef STATISTICS_RECORDER_HPP
#define STATISTICS_RECORDER_HPP

#include <fstream>
#include <vector>

// Define the BucketStatistics structure
struct BucketStatistics
{
    int positive_cases = 0;    // Count of positive cases
    int negative_cases = 0;    // Count of negative cases
    int total_pairs = 0;       // Total number of read pairs
    int bucket_size = 0;       // Size of the bucket
};

// Define the StatisticsRecorder class
class StatisticsRecorder
{
public:
    std::vector<BucketStatistics> bucket_stats;

    // Constructor to initialize the bucket statistics with a specified number of buckets
    StatisticsRecorder(size_t num_buckets)
    {
        bucket_stats.resize(num_buckets);
    }

    // Method to record the statistics for each bucket
    void record_bucket(int bucket_id, int positive_cases, int negative_cases, int total_pairs, int bucket_size)
    {
        bucket_stats[bucket_id] = {bucket_size, total_pairs, positive_cases, negative_cases};
    }

    // Method to write the statistics to a file
    void write_statistics_to_file(const std::string &filename, char delimiter = '\t')
    {
        std::ofstream out_file(filename);
        out_file << "BucketSize" << delimiter
                 << "TotalPairs" << delimiter
                 << "PositiveCases" << delimiter
                 << "NegativeCases" << '\n';

        for (const auto &stat : bucket_stats)
        {
            out_file << stat.bucket_size << delimiter
                     << stat.total_pairs << delimiter
                     << stat.positive_cases << delimiter
                     << stat.negative_cases << '\n';
        }
        out_file.close();
    }
};

#endif // STATISTICS_RECORDER_HPP
