// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <atomic>  // for atomic_fetch_add, atomic
#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/config.hpp>       // for config
#include <hibf/layout/graph.hpp> // for node_data
#include <hibf/misc/timer.hpp>   // for concurrent, timer

namespace seqan::hibf::build
{

/*!\brief Contains information used for building.
 * \ingroup hibf_build
 */
struct build_data
{
    std::atomic<size_t> ibf_number{};

    seqan::hibf::config const & config;

    layout::graph ibf_graph{};

    std::vector<double> fpr_correction{};

    // Timers do not copy the stored duration upon copy construction/assignment
    mutable concurrent_timer index_allocation_timer{};
    mutable concurrent_timer user_bin_io_timer{};
    mutable concurrent_timer merge_kmers_timer{};
    mutable concurrent_timer fill_ibf_timer{};

    size_t request_ibf_idx()
    {
        return std::atomic_fetch_add(&ibf_number, 1u);
    }
};

} // namespace seqan::hibf::build
