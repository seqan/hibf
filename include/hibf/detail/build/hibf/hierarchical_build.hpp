// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/contrib/robin_hood.hpp>

#include <lemon/list_graph.h> /// Must be first include.

#include <seqan3/search/dream_index/detail/build/hibf/build_data.hpp>
#include <seqan3/search/dream_index/detail/build/hibf/construct_ibf.hpp>
#include <seqan3/search/dream_index/detail/build/hibf/insert_into_ibf.hpp>
#include <seqan3/search/dream_index/detail/build/hibf/loop_over_children.hpp>
#include <seqan3/search/dream_index/detail/build/hibf/update_user_bins.hpp>

namespace seqan3::hibf
{

template <seqan3::data_layout data_layout_mode, typename config_type>
size_t hierarchical_build(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                          lemon::ListDigraph::Node const & current_node,
                          build_data<data_layout_mode, config_type> & data,
                          bool is_root)
{
    auto & current_node_data = data.node_map[current_node];

    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    robin_hood::unordered_flat_set<size_t> kmers{};

    auto compute_kmers = [&data, &kmers] (chopper_pack_record const & record)
    {
        for (auto && hash_sequence : data.hibf_config.input[record.user_bin_index])
            for (auto hash : hash_sequence)
                kmers.insert(hash);
    };

    auto initialise_max_bin_kmers = [&] (robin_hood::unordered_flat_set<size_t> & kmers,
                                         std::vector<int64_t> & ibf_positions,
                                         std::vector<int64_t> & filename_indices,
                                         lemon::ListDigraph::Node const & node,
                                         build_data<data_layout_mode, config_type> & data)
        -> size_t
    {
        auto & node_data = data.node_map[node];

        if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
        {
            // recursively initialize favourite child first
            ibf_positions[node_data.max_bin_index] =
                hierarchical_build(kmers, node_data.favourite_child, data, false);
            return 1;
        }
        else // max bin is not a merged bin
        {
            // we assume that the max record is at the beginning of the list of remaining records.
            auto const & record = node_data.remaining_records[0];
            compute_kmers(record);
            update_user_bins(data, filename_indices, record);

            return record.number_of_bins.back();
        }
    };

    // initialize lower level IBF
    size_t const max_bin_tbs =
        initialise_max_bin_kmers(kmers, ibf_positions, filename_indices, current_node, data);

    auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, is_root);

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];

        if (is_root && record.number_of_bins.back() == 1) // no splitting needed
        {
            insert_into_ibf(data, record, ibf);
        }
        else
        {
            compute_kmers(record);
            insert_into_ibf(parent_kmers, kmers, record.number_of_bins.back(), record.bin_indices.back(), ibf, is_root);
        }

        update_user_bins(data, filename_indices, record);
        kmers.clear();
    }

    data.hibf->ibf_vector[ibf_pos] = std::move(ibf);
    data.hibf->next_ibf_id[ibf_pos] = std::move(ibf_positions);
    data.hibf->user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

    return ibf_pos;
}

} // namespace seqan3::hibf
