// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <atomic>
#include <hibf/std/new>

#include <hibf/config.hpp>
#include <hibf/detail/build/node_data.hpp>
#include <hibf/detail/timer.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

namespace hibf
{

struct build_data
{
    std::atomic<size_t> ibf_number{};

    config const & hibf_config;

    lemon::ListDigraph ibf_graph{};
    lemon::ListDigraph::NodeMap<node_data> node_map{ibf_graph};

    std::vector<double> fp_correction{};

    // Timers do not copy the stored duration upon copy construction/assignment
    mutable timer<concurrent::yes> wall_clock_timer{};
    mutable timer<concurrent::yes> bin_size_timer{};
    mutable timer<concurrent::yes> index_allocation_timer{};
    mutable timer<concurrent::yes> user_bin_io_timer{};
    mutable timer<concurrent::yes> merge_kmers_timer{};
    mutable timer<concurrent::yes> fill_ibf_timer{};
    mutable timer<concurrent::yes> store_index_timer{};

    size_t request_ibf_idx()
    {
        return std::atomic_fetch_add(&ibf_number, 1u);
    }
};

} // namespace hibf