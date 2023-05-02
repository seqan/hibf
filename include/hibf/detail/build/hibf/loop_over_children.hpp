// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <hibf/contrib/robin_hood.hpp>

#include <hibf/migration/execution_handler_parallel/execution_handler_parallel.hpp>

#include <lemon/list_graph.h> /// Must be first include.

#include <hibf/detail/build/hibf/build_data.hpp>
#include <hibf/detail/build/hibf/hierarchical_build.hpp>
#include <hibf/detail/build/hibf/insert_into_ibf.hpp>

namespace hibf
{

template <typename config_type>
void loop_over_children(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                        hibf::interleaved_bloom_filter<> & ibf,
                        std::vector<int64_t> & ibf_positions,
                        lemon::ListDigraph::Node const & current_node,
                        build_data<config_type> & data,
                        bool is_root)
{
    auto & current_node_data = data.node_map[current_node];
    std::vector<lemon::ListDigraph::Node> children{};

    for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
        children.emplace_back(data.ibf_graph.target(arc_it));

    if (children.empty())
        return;

    size_t const number_of_mutex = (data.node_map[current_node].number_of_technical_bins + 63) / 64;
    std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

    auto worker = [&](auto && index, auto &&)
    {
        auto & child = children[index];

        if (child != current_node_data.favourite_child)
        {
            robin_hood::unordered_flat_set<size_t> kmers{};
            size_t const ibf_pos = hierarchical_build(kmers, child, data, false);
            auto parent_bin_index = data.node_map[child].parent_bin_index;
            {
                size_t const mutex_id{parent_bin_index / 64};
                std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                ibf_positions[parent_bin_index] = ibf_pos;
                insert_into_ibf(kmers, 1, parent_bin_index, ibf);
                if (!is_root)
                    parent_kmers.insert(kmers.begin(), kmers.end());
            }
        }
    };

    size_t number_of_threads{};
    std::vector<size_t> indices(children.size());
    std::iota(indices.begin(), indices.end(), size_t{});

    if (is_root)
    {
        // Shuffle indices: More likely to not block each other. Optimal: Interleave
        std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
        number_of_threads = data.hibf_config.threads;
    }
    else
    {
        number_of_threads = 1u;
    }

    seqan3::detail::execution_handler_parallel executioner{number_of_threads};
    executioner.bulk_execute(std::move(worker), std::move(indices), []() {});
}

} // namespace hibf
