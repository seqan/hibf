// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cinttypes> // for uint64_t
#include <cstddef>   // for size_t

#include <hibf/build/build_data.hpp>         // for build_data
#include <hibf/contrib/robin_hood.hpp>       // for unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter
#include <hibf/layout/graph.hpp>

namespace seqan::hibf::build
{

/*!\brief Constructs an IBF of the HIBF.
 * \ingroup hibf_build
 */
seqan::hibf::interleaved_bloom_filter construct_ibf(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                                    robin_hood::unordered_flat_set<uint64_t> & kmers,
                                                    size_t const number_of_bins,
                                                    layout::graph::node const & node,
                                                    build_data & data,
                                                    bool is_root);

} // namespace seqan::hibf::build
