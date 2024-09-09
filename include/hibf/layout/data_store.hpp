// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>   // for assert
#include <cinttypes> // for uint64_t
#include <cstddef>   // for size_t
#include <numeric>   // for iota
#include <string>    // for string
#include <vector>    // for vector

#include <hibf/layout/layout.hpp>      // for layout
#include <hibf/misc/timer.hpp>         // for concurrent_timer
#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

namespace seqan::hibf::layout
{

/*!\brief Contains information used for the layout.
 * \ingroup hibf_layout
 */
struct data_store
{
    /*!\brief Stores information of the previous level of a given IBF.
     * \details
     * When computing a hierarchical layout, the `data_store` data structure is used with local copies for each IBF
     * within the hierarchical structure of the HIBF. To keep track of the hierarchy, the `previous_level` stores
     * information about the previous level (where the corresponding merged bin is located).
     */
    struct previous_level
    {
        std::vector<size_t> bin_indices{};
        std::string num_of_bins;

        bool empty() const
        {
            assert(bin_indices.empty() == num_of_bins.empty());
            return bin_indices.empty();
        }
    };

    /*!\name References to global instances of the HIBF.
     * \{
     */
    //!\brief The desired maximum false positive rate of the resulting index.
    double false_positive_rate{};

    //!\brief The layout that is built by layout::hierarchical_binning.
    layout * hibf_layout; // Will be modified by {simple,hierarchical}_binning.

    //!\brief The kmer counts associated with the above files used to layout user bin into technical bins.
    std::vector<size_t> const * kmer_counts{}; // Pointed to data should not be modified.

    //!\brief The hyperloglog sketches of all input files to estimate their size and similarities.
    std::vector<sketch::hyperloglog> const * sketches{}; // Pointed to data should not be modified.
    //!\}

    /*!\name Local Storage one IBF in the HIBF.
     * \details
     * These member variables change on each IBF of the HIBF s.t. the current IBF can be constructed from
     * the current subset of data. The same data is also used for the top level IBF that holds all the data.
     * \{
     */

    //!\brief The input is sorted and rearranged. To keep track without changing the input we store the positions.
    std::vector<size_t> positions = [this]()
    {
        std::vector<size_t> ps;
        ps.resize(this->kmer_counts->size());
        std::iota(ps.begin(), ps.end(), 0);
        return ps;
    }(); // GCOVR_EXCL_LINE

    //!\brief The false positive correction based on fp_rate, num_hash_functions and requested_max_tb.
    std::vector<double> fpr_correction{};

    //!\brief The correction factor for merged bins which are allowed to have a relaxed FPR.
    double relaxed_fpr_correction{};

    //!\brief Information about previous levels of the IBF if the algorithm is called recursively.
    previous_level previous{};

    //!\brief Matrix of estimates of merged bin cardinalites
    std::vector<uint64_t> union_estimates{};

    bool user_bins_arranged{false};

    //!\brief Tracks the time the algorithm spends on estimating the union of user bins (merged bins).
    concurrent_timer union_estimation_timer{};

    //!\brief Tracks the time the algorithm spends on rearranging user bins (merged bins).
    concurrent_timer rearrangement_timer{};
    //!\}
};

} // namespace seqan::hibf::layout
