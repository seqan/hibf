// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <algorithm>  // for __sort_fn, sort
#include <cassert>    // for assert
#include <cinttypes>  // for int64_t, uint16_t
#include <concepts>   // for unsigned_integral
#include <cstddef>    // for size_t
#include <functional> // for identity, less
#include <ranges>     // for forward_range, range, range_value_t
#include <utility>    // for addressof
#include <vector>     // for vector

#include <hibf/cereal/concepts.hpp>          // for cereal_archive
#include <hibf/config.hpp>                   // for config
#include <hibf/detail/timer.hpp>             // for concurrent, timer
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter
#include <hibf/layout/layout.hpp>            // for layout

#include <cereal/macros.hpp> // for CEREAL_SERIALIZE_FUNCTION_NAME

namespace seqan::hibf
{

/*!\brief The Hierarchical Interleaved Bloom Filter (HIBF) - Fast answers to set-membership queries for multiple bins.
 * \details
 *
 * This class improves the [seqan::hibf::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. It is especially suited if you want to index are many samples/sets/user bins
 * or if their sizes are unevenly distributed.
 *
 * Publication reference: https://doi.org/10.1186/s13059-023-02971-4
 *
 * ### Example
 *
 * \include test/snippet/hibf/hierarchical_interleaved_bloom_filter.cpp
 *
 * ### Cite
 *
 * *Mehringer, Svenja, et al. "Hierarchical Interleaved Bloom Filter: enabling ultrafast, approximate sequence queries."
 * Genome Biology 24.1 (2023): 1-25.*
 *
 * ## Constructing the HIBF
 *
 * The HIBF is constructed by passing a seqan::hibf::config. There are two options required to be set:
 * (1) seqan::hibf::config::input_fn and (2) seqan::hibf::config::number_of_user bins. For all other options
 * we have set sensible defaults.
 *
 * Here are all options with their defaults:
 *
 * \include test/snippet/hibf/hibf_construction.cpp
 *
 * Please see the documentation of seqan::hibf::config for details on how to configure the HIBF construction.
 *
 * ## Membership Queries with the HIBF
 *
 * To allow efficient, thread-safe membership queries, you need to use the
 * seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent.
 *
 * \include test/snippet/hibf/hierarchical_interleaved_bloom_filter.cpp
 *
 * You retrieve an membership_agent by calling seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent().
 *
 * You can then pass a **range of hashes** and a **threshold**.
 *
 * ### Thresholding
 *
 * Given a number `x` of hashes to query and a threshold value `t`, a query will return all user bin ids for which at
 * least `t` number of hashes have been found in the respective user bin in the HIBF index. In other words, the hit
 * count must be equal or greater than `t` (`count >= t`).
 *
 * For all practical applications it is recommended to research sensible thresholds based on the data, the false
 * positive rate, the length of the query and whether (canonical) k-mers, minimizers, syncmers,.. etc were used for
 * hashing genomic content.
 *
 * ## Counting Queries with the HIBF
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * seqan::hibf::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned seqan::hibf::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * # Details on the data structure
 *
 * The following gives some insights about the general design of the HIBF data structure. More details can be found
 * in the publication: https://doi.org/10.1186/s13059-023-02971-4
 *
 * ## Terminology
 *
 * ### Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs).
 *
 * ### User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * IBF, the user is interested in an answer that differentiates between these groups.
 *
 * ## Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [seqan::hibf::interleaved_bloom_filter][1], the user bins may be split across multiple technical
 * bins, or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores
 * another IBF that is built over the user bins constituting the merged bin. This lower-level IBF can then be used
 * to further distinguish between merged bins.
 *
 * In this example, user bin 1 was split into two technical bins. Bins 3, 4, and 5 were merged into a single technical
 * bin, and another IBF was added for the merged bin.
 * \image html hibf.svg width=90%
 *
 * The individual IBFs may have a different number of technical bins and differ in their sizes, allowing an efficient
 * distribution of the user bins.
 *
 * \see [seqan::hibf::interleaved_bloom_filter][1]
 */
class hierarchical_interleaved_bloom_filter
{
public:
    /*!\brief Manages membership queries for the seqan::hibf::hierarchical_interleaved_bloom_filter.
    * \see seqan::hibf::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
    * \details
    * In contrast to the [seqan::hibf::interleaved_bloom_filter][1], the result will consist of indices of user bins.
    */
    class membership_agent_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_bloom_filter() = default;                                              //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter const &) = default;                        //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    ~hierarchical_interleaved_bloom_filter() = default;            //!< Defaulted.

    /*!\brief Constructs the HIBF from the given configuration.
     * \details
     *
     * ## Configuration
     * Required options:
     * * `input_fn`
     * * `number_of_user_bins`
     *
     * Options recommended to adapt to your setup:
     * * `threads` - Choose number of threads depending on your hardware settings to speed up construction
     * * `maximum_false_positive_rate` - How many false positive answers can you tolerate? A low FPR (e.g. 0.001) is
     *   needed if you can tolerate a high RAM peak when using the HIBF but post-processing steps are heavy and FPs
     *   should be avoided. A high FPR (e.g. `0.3`) can be chosed if you want a very small HIBF and false positive
     *   can be easily filtered in the down-stream analysis
     *
     * ## Validation
     * \copybrief seqan::hibf::config::validate_and_set_defaults
     */
    hierarchical_interleaved_bloom_filter(config & configuration);

    /*!\brief [Advanced] Constructs the HIBF from a config and layout.
     * \details
     * This constructor makes it possible to avoid computing the layout on construction of an hibf by using a given
     * layout. A layout file can be constructed manually or via chopper (https://github.com/seqan/chopper)
     * or raptor-layout (https://github.com/seqan/raptor).
     */
    hierarchical_interleaved_bloom_filter(config & configuration, layout::layout const & layout);
    //!\}

    //!\brief The individual interleaved Bloom filters.
    std::vector<interleaved_bloom_filter> ibf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the ID of the lower
     * level IBF in ibf_vector.
     */
    std::vector<std::vector<int64_t>> next_ibf_id;

    /*!\brief Stores for each bin in each IBF of the HIBF the user bin ID.
    * \details
    * Assume we look up a bin `b` in IBF `i`, i.e. `ibf_bin_to_user_bin_id[i][b]`.
    * If `-1` is returned, bin `b` is a merged bin, there is no single user bin, we need to look into the
    * lower level IBF.
    * Otherwise, the returned value `j` is the corresponding user bin ID.
    */
    std::vector<std::vector<int64_t>> ibf_bin_to_user_bin_id{};

    //!\brief Returns a membership_agent to be used for counting.
    membership_agent_type membership_agent() const;

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan::hibf::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template <seqan::hibf::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf_vector);
        archive(next_ibf_id);
#ifdef RAPTOR_OLD_HIBF // Temporary compatibility with Raptor's HIBF.
        std::vector<std::string> filenames{};
        archive(filenames);
#endif
        archive(ibf_bin_to_user_bin_id);
    }

    /*!\name Timer
     * \brief Only contains values after the HIBF has been constructed.
     * \{
     */
    timer<concurrent::yes> index_allocation_timer{};
    timer<concurrent::yes> user_bin_io_timer{};
    timer<concurrent::yes> merge_kmers_timer{};
    timer<concurrent::yes> fill_ibf_timer{};
    //!\}
    //!\endcond
};

class hierarchical_interleaved_bloom_filter::membership_agent_type
{
private:
    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hierarchical_interleaved_bloom_filter const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive membership querying.
    template <std::ranges::forward_range value_range_t>
    void membership_for_impl(value_range_t && values, int64_t const ibf_idx, uint16_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(values);

        uint16_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];

            auto const current_filename_index = hibf_ptr->ibf_bin_to_user_bin_id[ibf_idx][bin];

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    membership_for_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                  // last bin
                     current_filename_index != hibf_ptr->ibf_bin_to_user_bin_id[ibf_idx][bin + 1]) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer.emplace_back(current_filename_index);
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent_type() = default;                                         //!< Defaulted.
    membership_agent_type(membership_agent_type const &) = default;            //!< Defaulted.
    membership_agent_type & operator=(membership_agent_type const &) = delete; //!< Deleted. hibf_ptr is const.
    membership_agent_type(membership_agent_type &&) = default;                 //!< Defaulted.
    membership_agent_type & operator=(membership_agent_type &&) = delete;      //!< Deleted. hibf_ptr is const.
    ~membership_agent_type() = default;                                        //!< Defaulted.

    /*!\brief Construct a membership_agent_type for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit membership_agent_type(hierarchical_interleaved_bloom_filter const & hibf) : hibf_ptr(std::addressof(hibf))
    {}
    //!\}

    //!\brief Stores the result of membership_for().
    std::vector<int64_t> result_buffer;

    //!\brief Sorts the results.
    void sort_results()
    {
        std::ranges::sort(result_buffer);
    }

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership for all user bins contained in this index, based on `values` and the `threshold`.
     * \param[in] values The values to process; must model std::ranges::forward_range.
     * \param[in] threshold Report a user bin if there are at least this many hits.
     * \returns A vector of user bin ids (index values) with successfull set membership query.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * Each value in `values` is queried against the index and all hits are accumulated. If the accumulated sum of hits
     * reaches the threshold for a user bin, that user bin (its index value) is returned.
     *
     * ### Example
     *
     * Lets assume that the hibf index is build on 3 user bins (UB_A, UB_B, and UB_C) and the user bins contain the
     * following hash values:
     *
     * * 0: UB_A = {4,5,6,11}
     * * 1: UB_B = {4,5,11,12}
     * * 2: UB_C = {4,5,6,7,9,10}
     *
     * Then the following query:
     * ```cpp
     *    auto agent = hibf.membership_agent();
     *    auto result = agent.membership_for(std::vector<size_t>{4,5,6,7}, 3); // result = {0,2}
     * ```
     * would return a vector that contains the index values 0 and 2, indicating that UB_A (hits 4,5,6) and
     * UB_C (hits 4,5,6,7) reached the threshold of `>= 3` hits. UB_B only counts 2 hits (hits 4,5) and is thus not
     * contained in the list of user bins with a successful query.
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & membership_for(value_range_t && values,
                                                              uint16_t const threshold) & noexcept
    {
        assert(hibf_ptr != nullptr);

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        membership_for_impl(values, 0, threshold);

        return result_buffer;
    }

    // `membership_for` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & membership_for(value_range_t && values,
                                                              uint16_t const threshold) && noexcept = delete;
    //!\}
};

inline hierarchical_interleaved_bloom_filter::membership_agent_type
hierarchical_interleaved_bloom_filter::membership_agent() const
{
    return typename hierarchical_interleaved_bloom_filter::membership_agent_type{*this};
}

} // namespace seqan::hibf
