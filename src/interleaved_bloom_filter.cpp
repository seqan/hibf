// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for clamp
#include <array>      // for array
#include <bit>        // for bit_ceil, countl_zero
#include <cstdint>    // for uint64_t
#include <cstring>    // for size_t, memcpy
#include <functional> // for equal_to, function
#include <stdexcept>  // for logic_error, invalid_argument
#include <vector>     // for vector

#include <hibf/build/bin_size_in_bits.hpp>   // for bin_size_in_bits
#include <hibf/config.hpp>                   // for config, insert_iterator
#include <hibf/contrib/robin_hood.hpp>       // for hash, unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_count, bin_index, bin_size, hash_...
#include <hibf/misc/bit_vector.hpp>          // for bit_vector
#include <hibf/misc/divide_and_ceil.hpp>     // for divide_and_ceil
#include <hibf/misc/unreachable.hpp>         // for assert, unreachable
#include <hibf/platform.hpp>                 // for HIBF_COMPILER_IS_GCC
#include <hibf/sketch/hyperloglog.hpp>       // for hyperloglog

namespace seqan::hibf
{

#if HIBF_COMPILER_IS_GCC
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wattributes"
#endif // HIBF_COMPILER_IS_GCC

interleaved_bloom_filter::interleaved_bloom_filter(seqan::hibf::bin_count const bins_,
                                                   seqan::hibf::bin_size const size,
                                                   seqan::hibf::hash_function_count const funs,
                                                   bool const track_occupancy_) :
    bins{bins_.value},
    bin_size_{size.value},
    hash_funs{funs.value},
    track_occupancy{track_occupancy_}
{
    if (bins == 0)
        throw std::logic_error{"The number of bins must be > 0."};
    if (hash_funs == 0 || hash_funs > 5)
        throw std::logic_error{"The number of hash functions must be > 0 and <= 5."};
    if (bin_size_ == 0)
        throw std::logic_error{"The size of a bin must be > 0."};

    hash_shift = std::countl_zero(bin_size_);
    bin_words = divide_and_ceil(bins, 64u);
    technical_bins = bin_words * 64u;
    resize(technical_bins * bin_size_);
    occupancy.resize(technical_bins, 0u);
}

size_t find_biggest_bin(config const & configuration)
{
    size_t bin_id{};
    size_t max_size{};
    seqan::hibf::sketch::hyperloglog sketch{configuration.sketch_bits};

#pragma omp parallel for schedule(dynamic) num_threads(configuration.threads) firstprivate(sketch)
    for (size_t i = 0u; i < configuration.number_of_user_bins; ++i)
    {
        sketch.reset();
        configuration.input_fn(i, insert_iterator{sketch});

        size_t const estimate = sketch.estimate();
#pragma omp critical
        {
            if (estimate > max_size)
            {
                max_size = estimate;
                bin_id = i;
            }
        }
    }

    return bin_id;
}

size_t max_bin_size(config & configuration, size_t const max_bin_elements)
{
    configuration.validate_and_set_defaults();

    size_t const max_size = [&]()
    {
        if (max_bin_elements != 0u)
            return max_bin_elements;

        // Use sketches to determine biggest bin.
        size_t const max_bin_id = find_biggest_bin(configuration);
        // Get exact count for biggest bin. Sketch estimate's accuracy depends on configuration.sketch_bits
        robin_hood::unordered_flat_set<uint64_t> kmers{};
        configuration.input_fn(max_bin_id, insert_iterator{kmers});
        return kmers.size();
    }();

    return build::bin_size_in_bits({.fpr = configuration.maximum_fpr, //
                                    .hash_count = configuration.number_of_hash_functions,
                                    .elements = max_size});
}

// config validation is done by max_bin_size
interleaved_bloom_filter::interleaved_bloom_filter(config & configuration, size_t const max_bin_elements) :
    interleaved_bloom_filter{seqan::hibf::bin_count{configuration.number_of_user_bins},
                             seqan::hibf::bin_size{max_bin_size(configuration, max_bin_elements)},
                             seqan::hibf::hash_function_count{configuration.number_of_hash_functions},
                             configuration.track_occupancy}
{
    size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(bin_count() / configuration.threads), 8u, 64u);

#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(configuration.threads)
    for (size_t i = 0u; i < configuration.number_of_user_bins; ++i)
    {
        configuration.input_fn(i, insert_iterator{*this, i});
    }
}

[[gnu::always_inline]] void interleaved_bloom_filter::emplace(size_t const value, bin_index const bin) noexcept
{
    assert(bin.value < bins);

    bool exists{track_occupancy};

    for (size_t i = 0; i < hash_funs; ++i)
    {
        size_t const idx = hash_and_fit(value, hash_seeds[i]) + bin.value;
        assert(idx < size());

        // Constructing the reference twice for tracking occupancy would impact performance.
        // No difference for emplace.
        seqan::hibf::bit_vector::reference bit_reference{(*this)[idx]};
        if (track_occupancy)
            exists &= bit_reference;
        bit_reference = true;
    };

    if (track_occupancy && !exists)
        ++occupancy[bin.value];
}

void interleaved_bloom_filter::clear(bin_index const bin) noexcept
{
    assert(bin.value < bins);
    for (size_t idx = bin.value, i = 0; i < bin_size_; idx += technical_bins, ++i)
        (*this)[idx] = 0;
}

bool interleaved_bloom_filter::try_increase_bin_number_to(seqan::hibf::bin_count const new_bin_count) noexcept
{
    size_t const new_bins = new_bin_count.value;
    size_t const new_bin_words = divide_and_ceil(new_bins, 64u);

    if (new_bins < bins || new_bin_words > bin_words)
        return false;

    bins = new_bins;
    return true;
}

void interleaved_bloom_filter::increase_bin_number_to(seqan::hibf::bin_count const new_bin_count)
{
    if (new_bin_count.value < bins)
        throw std::invalid_argument{"The number of new bins must be >= the current number of bins."};

    if (try_increase_bin_number_to(new_bin_count))
        return;

    size_t const new_bins = new_bin_count.value;
    size_t const new_bin_words = divide_and_ceil(new_bins, 64u);

    size_t const new_technical_bins = new_bin_words * 64u;
    size_t const new_bit_size = bin_size_ * new_technical_bins;
    size_t const old_bit_size = size();
    size_t const delta = new_technical_bins - technical_bins + 64;

    resize(new_bit_size);
    uint64_t * const ptr = data();

    //               old       new
    // |-------------|---------|
    // Backwards copy blocks of size (old_)technical_bins such that the new blocks are of size new_technical_bins.
    for (size_t new_block_end = new_bit_size, old_block_end = old_bit_size; old_block_end > 0;
         new_block_end -= new_technical_bins, old_block_end -= technical_bins)
    {
        size_t const stop = new_block_end - new_technical_bins;

        // Need to copy word-wise (64 bits) inside the block
        for (size_t i = new_block_end - delta, j = old_block_end - 64u; stop && i >= stop; i -= 64u, j -= 64u)
        {
            // We are working on the bit size, so we need to convert the indices to 64-bit pointer indices
            ptr[i / 64u] = ptr[j / 64u];
            ptr[j / 64u] = 0ULL;
        }
    }

    bins = new_bins;
    bin_words = new_bin_words;
    technical_bins = new_technical_bins;
    occupancy.resize(technical_bins, 0u);
}

[[gnu::always_inline]] bit_vector const &
interleaved_bloom_filter::containment_agent_type::bulk_contains(size_t const value) & noexcept
{
    assert(ibf_ptr != nullptr);
    if (result_buffer.size() != ibf_ptr->bin_count())
        seqan::hibf::unreachable();

    // Needed for auto-vectorization of loop. ibf_ptr->bin_words could change bewtween loops.
    size_t const bin_words_ = ibf_ptr->bin_words;
    size_t const hash_funs_ = ibf_ptr->hash_funs;

    // Removes case for bin_words_ == 0u. The same statment inside the switch-case wouldn't have that effect.
    if (bin_words_ == 0u)
        seqan::hibf::unreachable();
    if (hash_funs_ == 0u)
        seqan::hibf::unreachable();

    for (size_t i = 0; i < hash_funs_; ++i)
        bloom_filter_indices[i] = ibf_ptr->hash_and_fit(value, ibf_ptr->hash_seeds[i]) / 64u;

    uint64_t * const raw = result_buffer.data();
    uint64_t const * const ibf_data = ibf_ptr->data();
    std::memcpy(raw, ibf_data + bloom_filter_indices[0], sizeof(uint64_t) * bin_words_);

    // GCOVR_EXCL_START
    // NOLINTBEGIN(bugprone-branch-clone)
    auto impl = [&]<size_t extent = 0u>()
    {
        for (size_t i = 1; i < hash_funs_; ++i)
        {
            uint64_t const * const ibf_raw = ibf_data + bloom_filter_indices[i];

            if constexpr (extent == 0u)
            {
#pragma omp simd
                for (size_t j = 0; j < bin_words_; ++j)
                    raw[j] &= ibf_raw[j];
            }
            else if constexpr (extent == 2u || extent == 4u || extent == 8u)
            {
#pragma omp simd
                for (size_t j = 0; j < extent; ++j)
                    raw[j] &= ibf_raw[j];
            }
            else
            {
                for (size_t j = 0; j < extent; ++j)
                    raw[j] &= ibf_raw[j];
            }
        }
    };
    // NOLINTEND(bugprone-branch-clone)

    // https://godbolt.org/z/rqaeWGGer
    // Having the loop inside impl instead of around the switch/case is faster.
    switch (bin_words_)
    {
    case 1u: // 1 AND (64 bit)
        impl.operator()<1u>();
        break;
    case 2u: // 1 SSE4 instruction (128 bit)
        impl.operator()<2u>();
        break;
    case 3u: // 3 AND (64 bit) + Loop Unroll
        impl.operator()<3u>();
        break;
    case 4u: // 1 AVX2 instruction (256 bit)
        impl.operator()<4u>();
        break;
    case 5u: // 5 AND (64 bit) + Loop Unroll
        impl.operator()<5u>();
        break;
    case 6u: // 6 AND (64 bit) + Loop Unroll
        impl.operator()<6u>();
        break;
    case 7u: // 7 AND (64 bit) + Loop Unroll
        impl.operator()<7u>();
        break;
    case 8u: // 1 AVX512 instruction (512 bit)
        impl.operator()<8u>();
        break;
    default: // Auto vectorize. Might create different versions.
        impl();
    }
    // GCOVR_EXCL_STOP

    return result_buffer;
}

#if HIBF_COMPILER_IS_GCC
#    pragma GCC diagnostic pop
#endif // HIBF_COMPILER_IS_GCC

} // namespace seqan::hibf
