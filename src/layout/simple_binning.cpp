// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm> // for max, copy
#include <cassert>   // for assert
#include <cstddef>   // for size_t
#include <limits>    // for numeric_limits
#include <vector>    // for vector

#include <hibf/layout/data_store.hpp>     // for data_store
#include <hibf/layout/layout.hpp>         // for layout
#include <hibf/layout/simple_binning.hpp> // for simple_binning
#include <hibf/misc/divide_and_ceil.hpp>  // for divide_and_ceil

namespace seqan::hibf::layout
{

size_t simple_binning::execute()
{
    assert(data != nullptr);
    assert(num_technical_bins > 0u);
    assert(num_user_bins > 0u);

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
        size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * data->fpr_correction[i + 1]);
        matrix[i][0] = divide_and_ceil(corrected_ub_cardinality, i + 1u);
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
                    static_cast<size_t>(ub_cardinality * data->fpr_correction[(i - i_prime)]);
                size_t score =
                    std::max<size_t>(divide_and_ceil(corrected_ub_cardinality, i - i_prime), matrix[i_prime][j - 1]);

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
        size_t const number_of_bins = (trace_i - next_i);
        size_t const cardinality = (*data->kmer_counts)[data->positions[trace_j]];
        size_t const corrected_cardinality = static_cast<size_t>(cardinality * data->fpr_correction[number_of_bins]);
        size_t const cardinality_per_bin = divide_and_ceil(corrected_cardinality, number_of_bins);

        data->hibf_layout->user_bins.emplace_back(data->previous.bin_indices,
                                                  bin_id,
                                                  number_of_bins,
                                                  data->positions[trace_j]);

        if (cardinality_per_bin > max_size)
        {
            max_id = bin_id;
            max_size = cardinality_per_bin;
        }

        bin_id += number_of_bins;

        trace_i = trace[trace_i][trace_j];
        --trace_j;
    }
    ++trace_i; // because we want the length not the index. Now trace_i == number_of_bins
    size_t const cardinality = (*data->kmer_counts)[data->positions[0]];
    size_t const corrected_cardinality = static_cast<size_t>(cardinality * data->fpr_correction[trace_i]);
    // NOLINTNEXTLINE(clang-analyzer-core.DivideZero)
    size_t const cardinality_per_bin = divide_and_ceil(corrected_cardinality, trace_i);

    data->hibf_layout->user_bins.emplace_back(data->previous.bin_indices, bin_id, trace_i, data->positions[0]);

    if (cardinality_per_bin > max_size)
        max_id = bin_id;

    return max_id;
}

} // namespace seqan::hibf::layout
