// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <functional>
#include <iterator>

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/platform.hpp>

namespace hibf
{

using insert_iterator = std::insert_iterator<robin_hood::unordered_flat_set<uint64_t>>;

struct config
{
    /*!\name General Configuration
     * \{
     */
    //!\brief A lambda how to hash your input. TODO: Detailed docu needed!
    std::function<void(size_t const, insert_iterator &&)> input_fn;

    //!\brief Number of user bins
    size_t number_of_user_bins{};

    //!\brief The number of hash functions for the IBFs.
    size_t number_of_hash_functions{2};

    //!\brief The desired false positive rate of the IBFs.
    double maximum_false_positive_rate{0.05};

    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t threads{1u};
    //!\}

    /*!\name Layout Configuration
     * \{
     */
    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits{12};

    //!\brief The maximum number of technical bins on each IBF in the HIBF.
    uint16_t tmax{64};

    /*\brief A scaling factor to influence the amount of merged bins produced by the algorithm.
     *
     * The higher alpha, the more weight is added artificially to the low level IBFs and thus the optimal
     * solution will contain less merged bins in the end because it costs more to merge bins.
     */
    double alpha{1.2};

    //!\brief The maximal cardinality ratio in the clustering intervals.
    double max_rearrangement_ratio{0.5};

    //!\brief Whether to estimate the union of kmer sets to possibly improve the binning or not.
    bool disable_estimate_union{true};

    //!\brief Whether to do a second sorting of the bins which takes into account similarity or not.
    bool disable_rearrangement{true};
    //!\}

    /*!\name Build Configuration
     * \{
     */
    // Related to k-mers
    bool disable_cutoffs{false};

    //!\brief If given, no layout algorithm is esxecuted but the layout from file is used for building.
    std::filesystem::path layout_file{};

    // Related to IBF
    // bool compressed{false};
    //!\}
};

} // namespace hibf
