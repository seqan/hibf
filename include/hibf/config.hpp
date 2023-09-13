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

/*!\brief The configuration used to build an (H)IBF
 *
 * The configuration can be used to construct an HIBF or IBF.
 *
 * When constructing an IBF, only the members `General Configuration` are considered, layout parameters from
 * the section `HIBF Layout Configuration` are ignored.
 *
 * \note If an option is marked [REQUIRED], an error will be thrown on (H)IBF
 *       construction if it is not set. An option is marked RECOMMENDED_TO_ADAPT if we give a sensible default
 *       but still recommend from experience that it is worth thinking about adjusting it to your data.
 *
 * ## The HIBF takes too long to construct?
 *
 * Check the documentation of the following options that influence the runtime:
 * * seqan::hibf::config::threads
 * * seqan::hibf::config::max_rearrangement_ratio
 * * seqan::hibf::config::disable_estimate_union
 * * seqan::hibf::config::disable_rearrangement
 *
 * ## The HIBF memory consumption is too high?
 *
 * Check the documentation of the following options that influence the memory consumption:
 * * seqan::hibf::config::threads
 * * seqan::hibf::config::number_of_hash_functions
 * * seqan::hibf::config::maximum_false_positive_rate
 *
 * ## General Configuration
 *
 * ### seqan::hibf::config::input_fn [REQUIRED]
 *
 * \copydetails seqan::hibf::config::input_fn
 *
 * ### seqan::hibf::config::number_of_user_bins [REQURIED]
 *
 * \copydetails seqan::hibf::config::number_of_user_bins
 *
 * ### seqan::hibf::config::number_of_hash_functions
 *
 * \copydetails seqan::hibf::config::number_of_hash_functions
 *
 * ### seqan::hibf::config::maximum_false_positive_rate [RECOMMENDED_TO_ADAPT]
 *
 * \copydetails seqan::hibf::config::maximum_false_positive_rate
 *
 * ### seqan::hibf::config::threads [RECOMMENDED_TO_ADAPT]
 *
 * \copydetails seqan::hibf::config::threads
 *
 * ## HIBF Layout Configuration
 *
 * ### seqan::hibf::config::sketch_bits
 *
 * \copydetails seqan::hibf::config::sketch_bits
 *
 * ### seqan::hibf::config::tmax
 *
 * \copydetails seqan::hibf::config::tmax
 *
 * ### seqan::hibf::config::max_rearrangement_ratio
 *
 * \copydetails seqan::hibf::config::max_rearrangement_ratio
 *
 * ### seqan::hibf::config::alpha
 *
 * \copydetails seqan::hibf::config::alpha
 *
 * ### seqan::hibf::config::disable_estimate_union
 *
 * \copydetails seqan::hibf::config::disable_estimate_union
 *
 * ### seqan::hibf::config::disable_rearrangement
 *
 * \copydetails seqan::hibf::config::disable_rearrangement
 *
 */
struct config
{
    /*!\name General Configuration
     * \{
     */
    /*!\brief A function for how to hash your input [REQUIRED].
     *
     * \attention This option is required!
     *
     * To efficiently construct the hierarchical structure of the HIBF, the input needs to be given as a
     * [function object](https://en.cppreference.com/w/cpp/utility/functional). The IBF construction can be done
     * in another way (see seqan::hibf::interleaved_bloom_filter) but has the config construction for consistency.
     *
     * \include test/snippet/hibf/config_input_fn_dummy.cpp
     *
     * It is important that the function object, a lambda in the example, has exactly two parameters:
     * (1) A `size_t const` which reflects the user bin ID to add hashes/values to in the (H)IBF.
     * (2) A `seqan::hibf::insert_iterator` that inserts hashes via `it = hash` where `hash` must be convertible to
     *     `uint64_t`.
     *
     * The above example would just insert a single hash (`42`) into each user bin.
     *
     * Let's look at two more realistic examples.
     *
     * **Inserting from a vector of values**:
     *
     * \include test/snippet/hibf/config_input_fn_vector.cpp
     *
     * **Inserting from a file**:
     *
     * \snippet test/snippet/hibf/config_input_fn_file.cpp main
     *
     * \attention Since the number of user bins cannot be inferred from the function object, you need to pass them by a
     * separate config member: `seqan::hibf::config::number_of_user_bins`.
     *
     */
    std::function<void(size_t const, insert_iterator &&)> input_fn{};

    /*!\brief The number of user bins.
     *
     * Since the data to construct the (H)IBF is given by a function object `seqan::hibf::config::input_fn`,
     * the number of user bins to consider must be given via this option.
     *
     * \include test/snippet/hibf/config_number_of_user_bins.cpp
     *
     * In this example, `12` user bins would be inserted into the (H)IBF, each only storing the hash `42`.
     */
    size_t number_of_user_bins{};

    /*!\brief The number of hash functions for the underlying Bloom Filters.
     *
     * The (H)IBF is based on the [Bloom Filter](https://en.wikipedia.org/wiki/Bloom_filter) data structure which
     * requires a set of hash functions. The number of hash functions used influences the speed and space but the
     * optimal number of hash functions is data dependent (see [Bloom Filter Calculator](https://hur.st/bloomfilter/)).
     *
     * Based on our experiments, we recommend a value of 2 (default).
     *
     * Be sure to experiment with this option with your data before changing it.
     */
    size_t number_of_hash_functions{2};

    /*!\brief The desired maximum false positive rate of the underlying Bloom Filters. [RECOMMENDED_TO_ADAPT]
     *
     * We ensure that when querying a single hash value in the (H)IBF, the probability of getting a false positive answer
     * will not exceed the value set for seqan::hibf::config::maximum_false_positive_rate.
     * The internal Bloom Filters will be configured accordingly. Individual Bloom Filters might have a different
     * but always lower false positive rate (FPR).
     *
     * Value must be in range [0,1].
     * Recommendation: default value (0.05)
     *
     * The FPR influences the memory consumption of the (H)IBF:
     * * A lower FPR limits the number of false-positive results, but increases the index size.
     * * A higher FPR can help to reduce memory consumption in cases where false-positive answers have little effect.
     *
     * \sa [Bloom Filter Calculator](https://hur.st/bloomfilter/).
     */
    double maximum_false_positive_rate{0.05};

    /*!\brief The number of threads to use during construction. [RECOMMENDED_TO_ADAPT]
     *
     * Using more threads increases the memory consumption during construction because the threads hold local
     * data. It can be beneficial to try a lower number of threads if you have limited RAM but many threads.
     *
     * Currently, the following parts of the HIBF construction process are parallelized:
     * * Sketch computation during layouting. Each thread processes a single user bin, keeping the total amount of
     *   hashes in memory before computing the sketch.
     * * The hierarchical, bottom-up build process. Lower level IBFs of the HIBF that are independent of each other
     *   may be built in parallel. The threads hold local hash sets needed for higher level merged bins.
     */
    size_t threads{1u};
    //!\}

    /*!\name HIBF Layout Configuration
     * \{
     */
    /*!\brief The number of bits for HyperLogLog sketches.
     *
     * The HIBF layout algorithm estimates the user bin sizes by computing
     * [HyperLogLog](https://en.wikipedia.org/wiki/HyperLogLog) sketches.
     *
     * A key parameter of HyperLogLog sketches is the number of bits of a hash value to consider for sketching.
     * Fewer bits accelerate the sketching process but decrease the accuracy.
     *
     * Value must be in range [5,32].
     * Recommendation: default value (12).
     *
     * Be sure to experiment with this option with your data before changing it.
     */
    uint8_t sketch_bits{12};

    /*!\brief The maximum number of technical bins of each IBF in the HIBF.
     *
     * One of the key methods of the HIBF for increasing query speed is limiting the number of (technical) bins of the
     * IBF data structure used within the HIBF.
     *
     * Choosing a good tmax is not trivial.
     *
     * The smaller tmax, the more levels the layout needs to represent the data. This results in a higher space
     * consumption of the index. While querying each individual level is cheap, querying many levels might also lead to
     * an increased runtime.
     *
     * A good tmax is usually the square root of the number of user bins/samples rounded to the next multiple of 64.
     * Note that your tmax will always be rounded to the next multiple of 64.
     *
     * Value must be in range [0,max(size_t)].
     * Recommendation: default value (the default will compute ≈sqrt(samples)).
     */
    size_t tmax{};

    /*\brief A scaling factor to influence the amount of merged bins produced by the layout algorithm.
     *
     * The layout algorithm optimizes the space consumption of the resulting HIBF, but currently has no means of
     * optimizing the runtime for querying such an HIBF. In general, the ratio of merged bins and split bins influences
     * the query time because a merged bin always triggers another search on a lower level. To influence this ratio,
     * alpha can be used.
     *
     * The higher alpha, the less merged bins are chosen in the layout. This improves query times but leads to a bigger
     * index.
     *
     * Value must be in range [0,max(double)].
     * Recommendation: default value (1.2).
     * disable_estimate_union
     * Be sure to experiment with this option with your data before changing it.
     */
    double alpha{1.2};

    /*!\brief The maximal cardinality ratio in the clustering intervals of the layout rearrangement algorithm.
     *
     * This option can influence the layout rearrangement algorithm.  The layout rearrangement algorithm improves the
     * layout by rearranging user bins with similar content into close proximity s.t. their shared hash values
     * reduces the overall index size.
     * The algorithm only rearranges the order of user bins in fixed intervals. The higher –max-rearrangement-ratio,
     * the larger the intervals. This potentially improves the layout, but increases the runtime of the layout
     * algorithm.
     *
     * Value must be in range [0,1].
     * Recommendation: default value (0.5).
     */
    double max_rearrangement_ratio{0.5};

    /*!\brief Whether to disable union estimate of user bins to improve the layout.
     *
     * The layout algorithm usually estimates the union size of user bins based on their HyperLogLog sketches in order
     * to precisely estimate the resulting HIBF memory consumption. It improves the layout quality but is
     * computationally expensive.
     *
     * If you are constructing an HIBF with more than 100,000 samples, we recommend considering disabling this feature
     * for faster HIBF construction time.
     */
    bool disable_estimate_union{false};

    /*!\brief Whether to disable rearranging user bins based on their content similarity.
     *
     * The layout rearrangement algorithm improves the layout by rearranging user bins with similar content into close
     * proximity s.t. their shared hash values reduce the overall index size. It improves the layout quality but is
     * computationally expensive.
     *
     * If you are constructing an HIBF with more than 100,000 samples, we recommend considering disabling this feature
     * for faster HIBF construction time.
     */
    bool disable_rearrangement{false};
    //!\}

    void read_from(std::istream & stream);
    void write_to(std::ostream & stream) const;

    /*!\brief Checks several variables of seqan::hibf::config and sets default values if necessary.
     *
     * Required options:
     *   * seqan::hibf::config::number_of_user_bins must be set to a value other than `0`.
     *   * seqan::hibf::config::input_fn must be set.
     *
     * Constrains:
     *   * seqan::hibf::config::number_of_hash_functions must be in `[1,5]`.
     *   * seqan::hibf::config::maximum_false_positive_rate must be in `[0.0,1.0]`.
     *   * seqan::hibf::config::threads must be greater than `0`.
     *   * seqan::hibf::config::sketch_bits must be in `[5,32]`.
     *   * seqan::hibf::config::tmax must be at most `18446744073709551552`.
     *   * seqan::hibf::config::alpha must be positive.
     *   * seqan::hibf::config::max_rearrangement_ratio must be in `[0.0,1.0]`.
     *
     * Modifications:
     *   * Enabling seqan::hibf::config::disable_estimate_union or setting seqan::hibf::config::max_rearrangement_ratio
     *     to `0.0` also enables seqan::hibf::config::disable_rearrangement.
     *   * Not setting seqan::hibf::config::tmax, or setting it to `0`, results in a default tmax
     *     `std::ceil(std::sqrt(number_of_user_bins))` being used.
     *   * seqan::hibf::config::tmax is increased to the next multiple of 64.
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
