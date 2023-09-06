// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes>  // for uint16_t, uint32_t, uint64_t, uint8_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iosfwd>     // for istream, ostream
#include <iterator>   // for insert_iterator

#include <hibf/cereal/path.hpp>        // IWYU pragma: keep
#include <hibf/contrib/robin_hood.hpp> // for unordered_flat_set
#include <hibf/platform.hpp>

#include <cereal/access.hpp> // for access
#include <cereal/cereal.hpp> // for make_nvp, CEREAL_NVP

namespace seqan::hibf
{

using insert_iterator = std::insert_iterator<robin_hood::unordered_flat_set<uint64_t>>;

struct config
{
    /*!\name General Configuration
     * \{
     */
    //!\brief A lambda how to hash your input. TODO: Detailed docu needed!
    std::function<void(size_t const, insert_iterator &&)> input_fn{};

    //!\brief Number of user bins
    size_t number_of_user_bins{};

    //!\brief The number of hash functions for the IBFs.
    size_t number_of_hash_functions{2};

    //!\brief The desired false positive rate of the IBFs.
    double maximum_false_positive_rate{0.05};

    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t threads{1u};
    //!\}

    /*!\name HIBF Layout Configuration
     * \{
     */
    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits{12};

    //!\brief The maximum number of technical bins on each IBF in the HIBF.
    uint16_t tmax{};

    /*\brief A scaling factor to influence the amount of merged bins produced by the algorithm.
     *
     * The higher alpha, the more weight is added artificially to the low level IBFs and thus the optimal
     * solution will contain less merged bins in the end because it costs more to merge bins.
     */
    double alpha{1.2};

    //!\brief The maximal cardinality ratio in the clustering intervals.
    double max_rearrangement_ratio{0.5};

    //!\brief Whether to estimate the union of kmer sets to possibly improve the binning or not.
    bool disable_estimate_union{false};

    //!\brief Whether to do a second sorting of the bins which takes into account similarity or not.
    bool disable_rearrangement{false};
    //!\}

    void read_from(std::istream & stream);
    void write_to(std::ostream & stream) const;

    /*!\brief Checks several variables of seqan::hibf::config and sets default values if necessary.
     *
     * The following checks are performed and will throw an exception if the checks fail:
     * * seqan::hibf::config::number_of_user_bins must be greather than `0`.
     * * If seqan::hibf::config::tmax is `0`, seqan::hibf::config::number_of_user_bins must be
     *   smaller than `1ULL << 32`.
     *
     * The configuration might be modified before being passed to the HIBF construction algorithm:
     * * If seqan::hibf::config::disable_estimate_union is `true`, seqan::hibf::config::disable_rearrangement will be
     *   set to `true` . Without union estimation, no rearrangement can be done.
     * * If seqan::hibf::config::tmax is `0`, the default value of `std::ceil(std::sqrt(cfg.number_of_user_bins))` will
     *   be used.
     * * If seqan::hibf::config::tmax is **not** `0` but also not a multiple of 64, it is increased to the next multiple
     *   of 64. E.g., the value `60` will be increased to `64`, and `1000` to `1024`.
     */
    void validate_and_set_defaults();

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t version{1};
        archive(CEREAL_NVP(version));

        archive(CEREAL_NVP(number_of_user_bins));
        archive(CEREAL_NVP(number_of_hash_functions));
        archive(CEREAL_NVP(maximum_false_positive_rate));
        archive(CEREAL_NVP(threads));

        archive(CEREAL_NVP(sketch_bits));
        archive(CEREAL_NVP(tmax));
        archive(CEREAL_NVP(alpha));
        archive(CEREAL_NVP(max_rearrangement_ratio));
        archive(CEREAL_NVP(disable_estimate_union));
        archive(CEREAL_NVP(disable_rearrangement));
    }
};

} // namespace seqan::hibf
