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
    for (auto chunk : kmers | seqan::stl::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan::hibf::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
            ibf.emplace(value, bin_idx);
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
    for (auto && value : values)
        ibf.emplace(value, bin_index);
    local_fill_ibf_timer.stop();
    data.fill_ibf_timer += local_fill_ibf_timer;
}

} // namespace seqan::hibf::build
