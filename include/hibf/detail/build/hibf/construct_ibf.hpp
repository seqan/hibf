// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <lemon/list_graph.h> /// Must be first include.

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/build/hibf/bin_size_in_bits.hpp>
#include <hibf/detail/build/hibf/build_data.hpp>

namespace hibf
{

hibf::interleaved_bloom_filter construct_ibf(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                             robin_hood::unordered_flat_set<uint64_t> & kmers,
                                             size_t const number_of_bins,
                                             lemon::ListDigraph::Node const & node,
                                             build_data & data,
                                             bool is_root);

} // namespace hibf
