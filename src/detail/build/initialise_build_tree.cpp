// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements hibf::initialise_build_tree.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#include <lemon/core.h>       // for INVALID
#include <lemon/list_graph.h> // for ListDigraph

#include <utility> // for move

#include <hibf/detail/build/initialise_build_tree.hpp>    // for initialise_build_tree
#include <hibf/detail/build/node_data.hpp>                // for node_data
#include <hibf/detail/build/update_content_node_data.hpp> // for update_content_node_data
#include <hibf/detail/build/update_header_node_data.hpp>  // for update_header_node_data
#include <hibf/detail/layout/layout.hpp>                  // for layout

namespace hibf
{

void initialise_build_tree(layout::layout & hibf_layout,
                           lemon::ListDigraph & ibf_graph,
                           lemon::ListDigraph::NodeMap<node_data> & node_map)
{
    // Add high level node
    auto high_level_node = ibf_graph.addNode(); // high-level node = root node
    node_map.set(high_level_node, {0, hibf_layout.top_level_max_bin_id, 0, lemon::INVALID, {}});

    update_header_node_data(std::move(hibf_layout.max_bins), ibf_graph, node_map);
    update_content_node_data(std::move(hibf_layout.user_bins), ibf_graph, node_map);
}

} // namespace hibf
