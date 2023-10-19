// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm> // for max, copy
#include <cassert>   // for assert
#include <cstddef>   // for size_t
#include <optional>  // for nullopt, optional
#include <ranges>    // for __prev_fn, prev
#include <utility>   // for move
#include <vector>    // for vector

#include <hibf/layout/graph.hpp>  // for graph
#include <hibf/layout/layout.hpp> // for layout

namespace seqan::hibf::layout
{

void update_header_node_data(std::vector<layout::layout::max_bin> const & header_max_bins,
                             seqan::hibf::layout::graph & ibf_graph)
{
#ifndef NDEBUG
    // head_max_bins must be sorted ascending by the number of bin indices (corresponds to the IBF levels)
    layout::layout::max_bin current = (header_max_bins.empty()) ? layout::layout::max_bin{} : *header_max_bins.begin();
    for (auto it = header_max_bins.begin() + 1; it < header_max_bins.end(); ++it)
    {
        assert(current.previous_TB_indices.size() <= it->previous_TB_indices.size());
        current = *it;
    }
#endif

    for (auto const & [bin_indices, max_id] : header_max_bins)
    {
        // we assume that the header lines are in the correct order
        // go down the tree until you find the matching parent
        seqan::hibf::layout::graph::node * parent = &ibf_graph.root; // start at root

        assert(!bin_indices.empty());
        auto bin_indices_sentinel = std::ranges::prev(bin_indices.end());
        for (auto it = bin_indices.begin(); it != bin_indices_sentinel; ++it)
        {
            for (seqan::hibf::layout::graph::node & child : parent->children)
            {
                if (child.parent_bin_index == *it)
                {
                    parent = &child;
                    break;
                }
            }
        }

        seqan::hibf::layout::graph::node new_child{{}, bin_indices.back(), max_id, 0u, std::nullopt, {}};
        parent->children.push_back(new_child);

        if (parent->max_bin_index == bin_indices.back())
            parent->favourite_child_idx = parent->children.size() - 1; // set favourite child to current child
    }
}

void update_content_node_data(std::vector<layout::layout::user_bin> const & layout_user_bins,
                              seqan::hibf::layout::graph & ibf_graph)
{
    // parse lines
    // -------------------------------------------------------------------------
    for (size_t user_bin = 0; user_bin < layout_user_bins.size(); ++user_bin)
    {
        auto const & record = layout_user_bins[user_bin];

        // go down the tree until you find the matching parent
        seqan::hibf::layout::graph::node * current_node = &ibf_graph.root; // start at root

        for (size_t i = 0; i < record.previous_TB_indices.size(); ++i)
        {
            size_t const bin = record.previous_TB_indices[i];
            size_t const num_tbs = 1;

            // update number of technical bins in current_node-IBF
            current_node->number_of_technical_bins = std::max(current_node->number_of_technical_bins, bin + num_tbs);

#ifndef NDEBUG
            bool found_next_node{false}; // sanity check
#endif
            for (auto & child : current_node->children)
            {
                if (child.parent_bin_index == bin)
                {
                    current_node = &child;
#ifndef NDEBUG
                    found_next_node = true;
#endif
                    break;
                }
            }
            assert(found_next_node);
        }

        size_t const bin = record.storage_TB_id;
        size_t const num_tbs = record.number_of_technical_bins;

        // update number of technical bins in current_node-IBF
        current_node->number_of_technical_bins = std::max(current_node->number_of_technical_bins, bin + num_tbs);

        if (record.storage_TB_id == current_node->max_bin_index)
            current_node->remaining_records.insert(current_node->remaining_records.begin(), record);
        else
            current_node->remaining_records.emplace_back(record);
    }
}

graph::graph(seqan::hibf::layout::layout const & hibf_layout)
{
    // initialise root node
    root.children = {};
    root.parent_bin_index = 0;
    root.max_bin_index = hibf_layout.top_level_max_bin_id;
    root.number_of_technical_bins = 0;
    root.favourite_child_idx = std::nullopt; // not known yet
    root.remaining_records = {};

    update_header_node_data(hibf_layout.max_bins, *this);
    update_content_node_data(hibf_layout.user_bins, *this);
}

} // namespace seqan::hibf::layout
