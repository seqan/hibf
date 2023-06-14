#include <algorithm> // for max
#include <limits>    // for numeric_limits

#include <hibf/detail/layout/layout.hpp> // for layout
#include <hibf/detail/layout/simple_binning.hpp>

namespace hibf::layout
{

size_t simple_binning::execute()
{
    assert(data != nullptr);

    std::vector<std::vector<size_t>> matrix(num_technical_bins); // rows
    for (auto & v : matrix)
        v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

    std::vector<std::vector<size_t>> trace(num_technical_bins); // rows
    for (auto & v : trace)
        v.resize(num_user_bins, std::numeric_limits<size_t>::max()); // columns

    size_t const extra_bins = num_technical_bins - num_user_bins + 1;

    // initialize first column (first row is initialized with inf)
    double const ub_cardinality = static_cast<double>((*data->kmer_counts)[data->positions[0]]);
    for (size_t i = 0; i < extra_bins; ++i)
    {
        size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * data->fp_correction[i + 1]);
        matrix[i][0] = corrected_ub_cardinality / (i + 1);
    }

    // we must iterate column wise
    for (size_t j = 1; j < num_user_bins; ++j)
    {
        double const ub_cardinality = static_cast<double>((*data->kmer_counts)[data->positions[j]]);

        for (size_t i = j; i < j + extra_bins; ++i)
        {
            size_t minimum{std::numeric_limits<size_t>::max()};

            for (size_t i_prime = j - 1; i_prime < i; ++i_prime)
            {
                size_t const corrected_ub_cardinality =
                    static_cast<size_t>(ub_cardinality * data->fp_correction[(i - i_prime)]);
                size_t score = std::max<size_t>(corrected_ub_cardinality / (i - i_prime), matrix[i_prime][j - 1]);

                // std::cout << "j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                minimum = (score < minimum) ? (trace[i][j] = i_prime, score) : minimum;
            }

            matrix[i][j] = minimum;
        }
    }

    // print_matrix(matrix, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());
    // print_matrix(trace, num_technical_bins, num_user_bins, std::numeric_limits<size_t>::max());

    // backtracking
    size_t trace_i = num_technical_bins - 1;
    size_t trace_j = num_user_bins - 1;

    size_t max_id{};
    size_t max_size{};

    size_t bin_id{};

    while (trace_j > 0)
    {
        size_t next_i = trace[trace_i][trace_j];
        size_t const kmer_count = (*data->kmer_counts)[data->positions[trace_j]];
        size_t const number_of_bins = (trace_i - next_i);
        size_t const kmer_count_per_bin = (kmer_count + number_of_bins - 1) / number_of_bins; // round up

        data->hibf_layout->user_bins.emplace_back(data->positions[trace_j],
                                                  data->previous.bin_indices,
                                                  number_of_bins,
                                                  bin_id);

        if (kmer_count_per_bin > max_size)
        {
            max_id = bin_id;
            max_size = kmer_count_per_bin;
        }

        bin_id += number_of_bins;

        trace_i = trace[trace_i][trace_j];
        --trace_j;
    }
    ++trace_i; // because we want the length not the index. Now trace_i == number_of_bins
    size_t const kmer_count = (*data->kmer_counts)[data->positions[0]];
    size_t const kmer_count_per_bin = (kmer_count + trace_i - 1) / trace_i;

    data->hibf_layout->user_bins.emplace_back(data->positions[0], data->previous.bin_indices, trace_i, bin_id);

    if (kmer_count_per_bin > max_size)
        max_id = bin_id;

    return max_id;
}

} // namespace hibf::layout
