#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <seqan3/search/dream_index/detail/configuration.hpp>
#include <seqan3/search/dream_index/detail/helper.hpp>
#include <seqan3/search/dream_index/detail/layout/hibf_statistics.hpp>
#include <seqan3/search/dream_index/detail/sketch/user_bin_sequence.hpp>

namespace seqan3::hibf
{

//!\brief Information about the previous IBF level to be passed down to ensure correct output.
struct previous_level
{
    std::vector<size_t> bin_index_values;
    std::string bin_indices;
    std::string num_of_bins;
    std::string estimated_sizes;
    std::string optimal_score;
    std::string correction;
    std::string tmax;

    bool empty() const
    {
        assert(((((bin_indices.empty() == num_of_bins.empty()) == estimated_sizes.empty()) == optimal_score.empty())
                == correction.empty())
               == tmax.empty());
        return bin_indices.empty();
    }
};

struct data_store
{
    /*!\name References to global instances of the HIBF.
     * \{
     */
    //!\brief The desired maximum false positive rate of the resulting index.
    double const false_positive_rate{};

    //!\brief A reference to the output stream to cache the results to.
    std::stringstream * output_buffer{nullptr};

    //!\brief A reference to the stream to cache the header to.
    std::stringstream * header_buffer{nullptr};

    //!\brief Stores what was formerly stored in the layout file header: For each merged bin, the max bin id.
    std::vector<std::pair<std::vector<size_t>, size_t>> * merged_bin_max_ids{};

    //!\brief An object starting at the top level IBF to collecting statistics about the HIBF on the way.
    hibf::hibf_statistics::level * stats{nullptr};
    //!\}

    /*!\name Local Storage one IBF in the HIBF.
     *
     * These member variables change on each IBF of the HIBF s.t. the current IBF can be constructed from
     * the current subset of data. The same data is also used for the top level IBF that holds all the data.
     * \{
     */
    //!\brief The file names of the user input. Since the input might be sorted, we need to keep track of the names.
    std::vector<std::string> filenames{};

    //!\brief The kmer counts associated with the above files used to layout user bin into technical bins.
    std::vector<size_t> kmer_counts{};

    //!\brief The hyperloglog sketches of all input files to estimate their size and similarities.
    std::vector<sketch::hyperloglog> sketches{};

    //!\brief Extra information in given in the input file
    std::vector<std::string> extra_information_strings{};

    //!\brief Extra information in given in the input file
    std::vector<std::vector<std::string>> extra_information{};

    //!\brief The false positive correction based on fp_rate, num_hash_functions and requested_max_tb.
    std::vector<double> fp_correction{};

    //!\brief Stores sketches if needed and provides utility functions for user bin rearrangement or union estimation.
    sketch::user_bin_sequence sketch_toolbox{};

    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    hibf::previous_level previous{};

    //!\brief Matrix of estimates of merged bin cardinalites
    std::vector<uint64_t> union_estimates{};

    bool user_bins_arranged{false};
    //!\}

    //!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
    void compute_fp_correction(double const fp_rate, size_t const num_hash_functions, size_t const requested_max_tb)
    {
        size_t const max_tb = next_multiple_of_64(requested_max_tb);

        fp_correction.resize(max_tb + 1, 0.0);
        fp_correction[1] = 1.0;

        double const denominator = std::log(1 - std::exp(std::log(fp_rate) / num_hash_functions));

        for (size_t i = 2; i <= max_tb; ++i)
        {
            double const tmp = 1.0 - std::pow(1 - fp_rate, static_cast<double>(i));
            fp_correction[i] = std::log(1 - std::exp(std::log(tmp) / num_hash_functions)) / denominator;
            assert(fp_correction[i] >= 1.0);
        }
    }
};

} // namespace seqan3
