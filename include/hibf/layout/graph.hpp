// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::layout::graph.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm> // for copy
#include <cstddef>   // for size_t
#include <optional>  // for nullopt, optional
#include <vector>    // for vector

#include <hibf/layout/layout.hpp> // for layout

namespace seqan::hibf::layout
{

/*!\brief Contains the layout graph structure.
 * \ingroup hibf_layout
 */
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
        std::optional<size_t> favourite_child_idx{std::nullopt};
        std::vector<layout::layout::user_bin> remaining_records{}; // non-merged bins (either split or single)

        bool max_bin_is_merged() const
        {
            return favourite_child_idx.has_value();
        }

        // Doesn't work, because the type is incomplete. To compare node, a comparison for the children member is needed.
        // But children is a std::vector<node>, so a comparison for node is needed to compare children.
        // https://godbolt.org/z/arrr4YKae
        // friend auto operator<=>(node const &, node const &) = default;
    };

    node root;
};

} // namespace seqan::hibf::layout
