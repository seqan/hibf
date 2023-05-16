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
#include <hibf/detail/build/hibf/insert_into_ibf.hpp>
#include <hibf/detail/build/hibf/update_parent_kmers.hpp>

namespace hibf
{

hibf::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                                 robin_hood::unordered_flat_set<uint64_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data & data,
                                                 bool is_root)
{
    auto & node_data = data.node_map[node];

    size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(kmers.size()) / number_of_bins))};
    double const bin_bits{static_cast<double>(bin_size_in_bits(kmers_per_bin, data.hibf_config.number_of_hash_functions, data.hibf_config.maximum_false_positive_rate))};
    hibf::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    hibf::bin_count const bin_count{node_data.number_of_technical_bins};

    timer<concurrent::no> local_index_allocation_timer{};
    local_index_allocation_timer.start();
    hibf::interleaved_bloom_filter<> ibf{bin_count, bin_size, hibf::hash_function_count{data.hibf_config.number_of_hash_functions}};
    local_index_allocation_timer.stop();
    data.index_allocation_timer += local_index_allocation_timer;

    insert_into_ibf(kmers, number_of_bins, node_data.max_bin_index, ibf, data.fill_ibf_timer);
    if (!is_root)
        update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);

    return ibf;
}

} // namespace hibf
