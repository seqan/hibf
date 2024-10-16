// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for equal_to, function
#include <iterator>   // for counted_iterator
#include <vector>     // for vector

#include <hibf/build/build_data.hpp>                // for build_data
#include <hibf/build/insert_into_ibf.hpp>           // for insert_into_ibf
#include <hibf/config.hpp>                          // for insert_iterator, config
#include <hibf/contrib/robin_hood.hpp>              // for Table, hash, unordered_flat_set
#include <hibf/contrib/std/chunk_view.hpp>          // for chunk, chunk_fn, chunk_view
#include <hibf/contrib/std/detail/adaptor_base.hpp> // for operator|
#include <hibf/interleaved_bloom_filter.hpp>        // for bin_index, interleaved_bloom_filter
#include <hibf/layout/layout.hpp>                   // for layout
#include <hibf/misc/divide_and_ceil.hpp>            // for divide_and_ceil
#include <hibf/misc/timer.hpp>                      // for serial_timer, concurrent_timer

namespace seqan::hibf::build
{

template <bool use_exists>
inline void
dispatch_emplace(seqan::hibf::interleaved_bloom_filter & ibf, auto && values, seqan::hibf::bin_index const bin_index)
{
    if constexpr (use_exists)
    {
        for (auto && value : values)
            ibf.emplace_exists(value, bin_index);
    }
    else
    {
        for (auto && value : values)
            ibf.emplace(value, bin_index);
    }
}

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(build_data const & data,
                     robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan::hibf::interleaved_bloom_filter & ibf,
                     concurrent_timer & fill_ibf_timer)
{
    // TODO:
    // 1) Some bins might have no element inserted.
    //    Example: Split 378 k-mers into 31 bins
    //             chunk_size will be 13.
    //             Bin sizes will be: 13,13,...,13,13,5,0
    // 2) It might not even be possible to just put 1 k-mer into each bin in extreme cases.
    //
    // We could get rid of the IBF's occupied_bins bit_vector, and just use the occupancy to determine empty bins if
    // it's guaranteed that every non-empty bin has a non-zero occupancy.
    size_t const chunk_size = divide_and_ceil(kmers.size(), number_of_bins);
    size_t chunk_number{};

    bool const use_exists = data.config.empty_bin_fraction > 0.0;

    serial_timer local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    auto chunk_view = seqan::stl::views::chunk(kmers, chunk_size);
    // For testing. Assert failes in `hibf_test.small_example_with_direct_hashes`.
#if 0
    assert(chunk_view.size() == number_of_bins
           || ((std::cerr << "chunk_view.size() = " << chunk_view.size() << "\nnumber_of_bins = " << number_of_bins
                          << '\n')
               && false));
#endif
    for (auto && chunk : chunk_view)
    {
        assert(chunk_number < number_of_bins);
        seqan::hibf::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        if (use_exists)
            dispatch_emplace<true>(ibf, std::move(chunk), bin_idx);
        else
            dispatch_emplace<false>(ibf, std::move(chunk), bin_idx);
    }
    local_fill_ibf_timer.stop();
    fill_ibf_timer += local_fill_ibf_timer;
}

void insert_into_ibf(build_data const & data,
                     layout::layout::user_bin const & record,
                     seqan::hibf::interleaved_bloom_filter & ibf)
{
    auto const bin_index = seqan::hibf::bin_index{static_cast<size_t>(record.storage_TB_id)};
    std::vector<uint64_t> values;

    serial_timer local_user_bin_io_timer{};
    local_user_bin_io_timer.start();
    data.config.input_fn(record.idx, insert_iterator{values});
    local_user_bin_io_timer.stop();
    data.user_bin_io_timer += local_user_bin_io_timer;

    serial_timer local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    if (data.config.empty_bin_fraction > 0.0)
        dispatch_emplace<true>(ibf, std::move(values), bin_index);
    else
        dispatch_emplace<false>(ibf, std::move(values), bin_index);
    local_fill_ibf_timer.stop();
    data.fill_ibf_timer += local_fill_ibf_timer;
}

} // namespace seqan::hibf::build
