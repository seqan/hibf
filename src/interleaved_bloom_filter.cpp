// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>  // for clamp, max
#include <array>      // for array
#include <bit>        // for bit_ceil, countl_zero
#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t
#include <cstring>    // for size_t, memcpy
#include <functional> // for function
#include <iterator>   // for inserter
#include <stdexcept>  // for logic_error, invalid_argument

#include <hibf/build/bin_size_in_bits.hpp>   // for bin_size_in_bits
#include <hibf/config.hpp>                   // for config
#include <hibf/contrib/robin_hood.hpp>       // for unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_count, bin_index, bin_size, hash_...
#include <hibf/platform.hpp>                 // for HIBF_COMPILER_IS_GCC

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

void interleaved_bloom_filter::emplace(size_t const value, bin_index const bin) noexcept
{
    assert(bin.value < bins);
    for (size_t i = 0; i < hash_funs; ++i)
    {
        size_t idx = hash_and_fit(value, hash_seeds[i]);
        idx += bin.value;
        assert(idx < data.size());
        data[idx] = 1;
    };
}

void interleaved_bloom_filter::clear(bin_index const bin) noexcept
{
    assert(bin.value < bins);
    for (size_t idx = bin.value, i = 0; i < bin_size_; idx += technical_bins, ++i)
        data[idx] = 0;
}

void interleaved_bloom_filter::increase_bin_number_to(seqan::hibf::bin_count const new_bins_)
{
    size_t new_bins = new_bins_.value;

    if (new_bins < bins)
        throw std::invalid_argument{"The number of new bins must be >= the current number of bins."};

    // Equivalent to ceil(new_bins / 64)
    size_t new_bin_words = (new_bins + 63) >> 6;

    bins = new_bins;

    if (new_bin_words == bin_words) // No need for internal resize if bin_words does not change.
        return;

    size_t new_technical_bins = new_bin_words << 6;
    size_t new_bits = bin_size_ * new_technical_bins;

    size_t idx_{new_bits}, idx{data.size()};
    size_t delta = new_technical_bins - technical_bins + 64;

    data.resize(new_bits);

    for (size_t i = idx_, j = idx; j > 0; i -= new_technical_bins, j -= technical_bins)
    {
        size_t stop = i - new_technical_bins;

        for (size_t ii = i - delta, jj = j - 64; stop && ii >= stop; ii -= 64, jj -= 64)
        {
            uint64_t old = data.get_int(jj);
            data.set_int(jj, 0);
            data.set_int(ii, old);
        }
    }

    bin_words = new_bin_words;
    technical_bins = new_technical_bins;
}

#if HIBF_COMPILER_IS_GCC
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wattributes"
#endif // HIBF_COMPILER_IS_GCC
[[gnu::always_inline]] binning_bitvector const &
interleaved_bloom_filter::membership_agent_type::bulk_contains(size_t const value) & noexcept
{
#if HIBF_COMPILER_IS_GCC
#    pragma GCC diagnostic pop
#endif // HIBF_COMPILER_IS_GCC
    assert(ibf_ptr != nullptr);
    assert(result_buffer.size() == ibf_ptr->bin_count());

    // Needed for auto-vectorization of loop. ibf_ptr->bin_words could change bewtween loops.
    size_t const bin_words = ibf_ptr->bin_words;
    size_t const hash_funs = ibf_ptr->hash_funs;

#ifndef NDEBUG
    assert(bin_words != 0u);
    assert(hash_funs != 0u);
#else
    // Removes case for bin_words == 0u. The same statment inside the switch-case wouldn't have that effect.
    if (bin_words == 0u)
        __builtin_unreachable();
    if (hash_funs == 0u)
        __builtin_unreachable();
#endif

    for (size_t i = 0; i < hash_funs; ++i)
        bloom_filter_indices[i] = ibf_ptr->hash_and_fit(value, ibf_ptr->hash_seeds[i]) >> 6;

    uint64_t * const raw = result_buffer.raw_data().data(); // TODO: std::assume_aligned<64> once memory-aligned
    uint64_t const * const ibf_data = ibf_ptr->data.data(); // TODO: std::assume_aligned<64> once memory-aligned
    std::memcpy(raw, ibf_data + bloom_filter_indices[0], sizeof(uint64_t) * bin_words);

    // https://godbolt.org/z/1nbhvqeGj
    // Having the loop inside is faster.
    // GCOVR_EXCL_START
    switch (bin_words)
    {
    case 1u: // 1 AND (64 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
            raw[0] &= ibf_raw[0];
        }
        break;
    case 2u: // 1 SSE4 instruction (128 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 2u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    case 3u: // 1 SSE4 instruction (128 bit) + 1 AND (64 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 3u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    case 4u: // 1 AVX2 instruction (256 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 4u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    case 5u: // 1 AVX2 instruction (256 bit) + 1 AND (64 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 5u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    case 6u: // 1 AVX2 instruction (256 bit) + 1 SSE4 instruction (128 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 6u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    case 7u: // 1 AVX2 instruction (256 bit) + 1 SSE4 instruction (128 bit) + 1 AND (64 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 7u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    case 8u: // 1 AVX512 instruction (512 bit)
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < 8u; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
        break;
    default: // Auto vectorize. Might create different versions.
        for (size_t i = 1; i < hash_funs; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];
#pragma omp simd
            for (size_t batch = 0; batch < bin_words; ++batch)
                raw[batch] &= ibf_raw[batch];
        }
    }
    // GCOVR_EXCL_STOP

    return result_buffer;
}

} // namespace seqan::hibf
