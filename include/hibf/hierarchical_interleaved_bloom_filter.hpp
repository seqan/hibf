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

#include <hibf/config.hpp>                   // for config
#include <hibf/detail/cereal/concepts.hpp>   // for cereal_archive
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter
#include <hibf/user_bins_type.hpp>           // for user_bins_type

#include <cereal/macros.hpp> // for CEREAL_SERIALIZE_FUNCTION_NAME

namespace hibf
{

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [hibf::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [hibf::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [hibf::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [hibf::interleaved_bloom_filter][1], the number of user bins
 * and technical bins is always the same. This causes performance degradation when there are many user bins or the user
 * bins are unevenly distributed.
 *
 * # Terminology
 *
 * ## Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs).
 *
 * ## User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * IBF, the user is interested in an answer that differentiates between these groups.
 *
 * # Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [hibf::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins
 * , or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores
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
 * ## Querying
 * To query the Hierarchical Interleaved Bloom Filter for values, call
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent() and use the returned
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent.
 * In contrast to the [hibf::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * hibf::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned hibf::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * [1]: https://docs.seqan.de/seqan/3.0.3/classseqan3_1_1interleaved__bloom__filter.html
 */
class hierarchical_interleaved_bloom_filter
{
public:
    /*!\brief Manages membership queries for the hibf::hierarchical_interleaved_bloom_filter.
    * \see hibf::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
    * \details
    * In contrast to the [hibf::interleaved_bloom_filter][1], the result will consist of indices of user bins.
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

    hierarchical_interleaved_bloom_filter(config const & configuration);
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

    //!\brief The underlying user bins.
    user_bins_type user_bins;

    //!\brief Returns a membership_agent to be used for counting.
    membership_agent_type membership_agent() const;

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy hibf::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     * \sa https://docs.seqan.de/seqan/3.2.0/group__io.html#serialisation
     */
    template <hibf::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf_vector);
        archive(next_ibf_id);
        archive(user_bins);
    }
    //!\endcond
};

class hierarchical_interleaved_bloom_filter::membership_agent_type
{
private:
    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hierarchical_interleaved_bloom_filter const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive membership querying.
    template <std::ranges::forward_range value_range_t>
    void bulk_contains_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(values);

        uint16_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];

            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_contains_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
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

    //!\brief Stores the result of bulk_contains().
    std::vector<int64_t> result_buffer;

    /*!\name Lookup
        * \{
        */
    /*!\brief Determines set membership of given values, and returns the user bin indices of occurrences.
        * \param[in] values The values to process; must model std::ranges::forward_range.
        * \param[in] threshold Report a user bin if there are at least this many hits.
        *
        * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
        * \attention Sequential calls to this function invalidate the previously returned reference.
        *
        * \details
        *
        * ### Thread safety
        *
        * Concurrent invocations of this function are not thread safe, please create a
        * hibf::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
        */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) & noexcept
    {
        assert(hibf_ptr != nullptr);

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        bulk_contains_impl(values, 0, threshold);

        std::ranges::sort(result_buffer); // TODO: necessary?

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values,
                                                             size_t const threshold) && noexcept = delete;
    //!\}
};

inline hierarchical_interleaved_bloom_filter::membership_agent_type
hierarchical_interleaved_bloom_filter::membership_agent() const
{
    return typename hierarchical_interleaved_bloom_filter::membership_agent_type{*this};
}

} // namespace hibf
