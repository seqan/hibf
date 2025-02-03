// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __fn, all_of, fill_n
#include <cassert>    // for assert
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t
#include <functional> // for equal_to, function
#include <iterator>   // for counted_iterator
#include <ranges>     // for subrange
#include <vector>     // for vector

#include <hibf/build/build_data.hpp>         // for build_data
#include <hibf/build/insert_into_ibf.hpp>    // for insert_into_ibf
#include <hibf/config.hpp>                   // for insert_iterator, config
#include <hibf/contrib/robin_hood.hpp>       // for Table, unordered_flat_set, hash
#include <hibf/contrib/std/chunk_view.hpp>   // for chunk_view, chunk, chunk_fn
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index
#include <hibf/layout/layout.hpp>            // for layout
#include <hibf/misc/divide_and_ceil.hpp>     // for divide_and_ceil
#include <hibf/misc/timer.hpp>               // for serial_timer, concurrent_timer

namespace seqan::hibf::build
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan::hibf::interleaved_bloom_filter & ibf,
                     concurrent_timer & fill_ibf_timer)
{
    size_t const chunk_size = divide_and_ceil(kmers.size(), number_of_bins);
    size_t chunk_number{};

    serial_timer local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    auto chunk_view = seqan::stl::views::chunk(kmers, chunk_size);
    for (auto && chunk : chunk_view)
    {
        assert(chunk_number < number_of_bins);
        seqan::hibf::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (auto && value : chunk)
            ibf.emplace(value, bin_idx);
    }

    assert(chunk_view.size() <= number_of_bins);
    // Edge case: If there are not enough k-mers to emplace at least one value into each bin, set the occupancy of
    // the left over bins to 1.
    // GCOVR_EXCL_START
    if (ibf.track_occupancy && chunk_view.size() < number_of_bins)
    {
        size_t const diff = number_of_bins - chunk_view.size();
        auto it = ibf.occupancy.begin() + bin_index + chunk_view.size();
        assert(std::ranges::all_of(it,
                                   it + diff,
                                   [](size_t value)
                                   {
                                       return value == 0u;
                                   }));
        std::ranges::fill_n(it, diff, 1u);
    }
    // GCOVR_EXCL_STOP

    local_fill_ibf_timer.stop();
    fill_ibf_timer += local_fill_ibf_timer;
}

void insert_into_ibf(build_data const & data,
                     layout::layout::user_bin const & record,
                     seqan::hibf::interleaved_bloom_filter & ibf)
{
    serial_timer local_user_bin_io_timer{};
    serial_timer local_fill_ibf_timer{};
    local_user_bin_io_timer.start();
    local_fill_ibf_timer.start();
    data.config.input_fn(record.idx, insert_iterator{ibf, record.storage_TB_id});
    local_user_bin_io_timer.stop();
    local_fill_ibf_timer.stop();
    data.user_bin_io_timer += local_user_bin_io_timer;
    data.fill_ibf_timer += local_fill_ibf_timer;
}

} // namespace seqan::hibf::build
