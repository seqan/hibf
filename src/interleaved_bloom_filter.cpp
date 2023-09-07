// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>  // for clamp, max
#include <bit>        // for bit_ceil, countl_zero
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iterator>   // for inserter
#include <stdexcept>  // for logic_error

#include <hibf/build/bin_size_in_bits.hpp>   // for bin_size_in_bits
#include <hibf/config.hpp>                   // for config
#include <hibf/contrib/robin_hood.hpp>       // for unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_count, bin_size, hash_function_count

#include <sdsl/int_vector.hpp> // for bit_vector

namespace seqan::hibf
{

interleaved_bloom_filter::interleaved_bloom_filter(seqan::hibf::bin_count bins_,
                                                   seqan::hibf::bin_size size,
                                                   seqan::hibf::hash_function_count funs)
{
    bins = bins_.value;
    bin_size_ = size.value;
    hash_funs = funs.value;

    if (bins == 0)
        throw std::logic_error{"The number of bins must be > 0."};
    if (hash_funs == 0 || hash_funs > 5)
        throw std::logic_error{"The number of hash functions must be > 0 and <= 5."};
    if (bin_size_ == 0)
        throw std::logic_error{"The size of a bin must be > 0."};

    hash_shift = std::countl_zero(bin_size_);
    bin_words = (bins + 63) >> 6;    // = ceil(bins/64)
    technical_bins = bin_words << 6; // = bin_words * 64
    data = sdsl::bit_vector(technical_bins * bin_size_);
}

size_t max_bin_size(config & configuration)
{
    configuration.validate_and_set_defaults();

    size_t max_size{};
    robin_hood::unordered_flat_set<uint64_t> kmers;
#pragma omp parallel for schedule(dynamic) num_threads(configuration.threads) private(kmers)
    for (size_t i = 0u; i < configuration.number_of_user_bins; ++i)
    {
        kmers.clear();
        configuration.input_fn(i, std::inserter(kmers, kmers.begin()));

#pragma omp critical
        max_size = std::max(max_size, kmers.size());
    }

    return build::bin_size_in_bits({.fpr = configuration.maximum_false_positive_rate,
                                    .hash_count = configuration.number_of_hash_functions,
                                    .elements = max_size});
}

// config validation is done by max_bin_size
interleaved_bloom_filter::interleaved_bloom_filter(config & configuration) :
    interleaved_bloom_filter{seqan::hibf::bin_count{configuration.number_of_user_bins},
                             seqan::hibf::bin_size{max_bin_size(configuration)},
                             seqan::hibf::hash_function_count{configuration.number_of_hash_functions}}
{
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(bin_count() / configuration.threads), 8u, 64u);
    robin_hood::unordered_flat_set<uint64_t> kmers;

#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(configuration.threads) private(kmers)
    for (size_t i = 0u; i < configuration.number_of_user_bins; ++i)
    {
        kmers.clear();
        configuration.input_fn(i, std::inserter(kmers, kmers.begin()));

        for (uint64_t const hash : kmers)
            emplace(hash, seqan::hibf::bin_index{i});
    }
}

} // namespace seqan::hibf
