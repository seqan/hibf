// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements seqan::hibf::layout::graph.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <optional>
#include <vector> // for vector, operator==

#include <hibf/layout/layout.hpp>

namespace seqan::hibf::layout
{

// Currently, the layout is structured by user bin.
struct graph
{

    /*!\name Constructors, destructor and assignment
     * \{
     */
    graph() = default;                          //!< Defaulted.
    graph(graph const &) = default;             //!< Defaulted.
    graph & operator=(graph const &) = default; //!< Defaulted.
    graph(graph &&) = default;                  //!< Defaulted.
    graph & operator=(graph &&) = default;      //!< Defaulted.
    ~graph() = default;                         //!< Defaulted.

    graph(seqan::hibf::layout::layout const & hibf_layout);
    //!\}

    struct node
    {
        std::vector<node> children{};

        size_t parent_bin_index{};
        size_t max_bin_index{};
        size_t number_of_technical_bins{};
        std::optional<uint32_t> favourite_child_idx{std::nullopt};
        std::vector<layout::layout::user_bin> remaining_records{}; // non-merged bins (either split or single)

        // Doesn't work, because the type is incomplete. To compare node, a comparison for the children member is needed.
        // But children is a std::vector<node>, so a comparison for node is needed to compare children.
        // https://godbolt.org/z/arrr4YKae
        // friend auto operator<=>(node const &, node const &) = default;
    };

    node root;
};

} // namespace seqan::hibf::layout
