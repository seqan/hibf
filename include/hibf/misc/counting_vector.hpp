// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::hibf::counting_vector.
 */

#pragma once

#include <algorithm>   // for fill
#include <array>       // for array
#include <bit>         // for countr_zero
#include <cassert>     // for assert
#include <cinttypes>   // for uint64_t, uint16_t
#include <concepts>    // for integral, same_as, unsigned_integral
#include <cstring>     // for size_t
#include <functional>  // for plus
#include <ranges>      // for range, forward_range, input_range, range_reference_t, range_value_t
#include <type_traits> // for remove_cvref_t
#include <utility>     // for addressof
#include <vector>      // for vector

#include <hibf/cereal/concepts.hpp>           // for cereal_archive
#include <hibf/config.hpp>                    // for config
#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/misc/bit_vector.hpp>           // for bit_vector

#include <cereal/macros.hpp>           // for CEREAL_SERIALIZE_FUNCTION_NAME
#include <cereal/types/base_class.hpp> // for base_class

namespace seqan::hibf
{

/*!\brief A data structure that behaves like a std::vector and can be used to consolidate the results of multiple calls
 *        to seqan::hibf::interleaved_bloom_filter::membership_agent_type::bulk_contains.
 * \ingroup ibf
 * \tparam value_t The type of the count. Must model std::integral.
 *
 * \details
 *
 * When using the seqan::hibf::interleaved_bloom_filter::membership_agent_type::bulk_contains operation, a common use
 * case is to add up, for example, the results for all k-mers in a query. This yields, for each bin, the number of
 * k-mers of a query that are in the respective bin. Such information can be used to apply further filtering or
 * abundance estimation based on the k-mer counts.
 *
 * The seqan::hibf::counting_vector offers an easy way to add up the individual
 * seqan::hibf::bit_vector by offering an `+=` operator.
 *
 * The `value_t` template parameter should be chosen in a way that no overflow occurs if all calls to `bulk_contains`
 * return a hit for a specific bin. For example, `uint8_t` will suffice when processing short Illumina reads, whereas
 * long reads will require at least `uint32_t`.
 *
 * ### Example
 *
 * \include test/snippet/ibf/counting_vector.cpp
 */
template <std::integral value_t>
class counting_vector : public std::vector<value_t>
{
private:
    //!\brief The base type.
    using base_t = std::vector<value_t>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_vector() = default;                                    //!< Defaulted.
    counting_vector(counting_vector const &) = default;             //!< Defaulted.
    counting_vector & operator=(counting_vector const &) = default; //!< Defaulted.
    counting_vector(counting_vector &&) = default;                  //!< Defaulted.
    counting_vector & operator=(counting_vector &&) = default;      //!< Defaulted.
    ~counting_vector() = default;                                   //!< Defaulted.

    using base_t::base_t;
    //!\}

    /*!\brief Bin-wise adds the bits of a seqan::hibf::bit_vector.
     * \param bit_vector The seqan::hibf::bit_vector.
     * \attention The counting_vector must be at least as big as `bit_vector`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/counting_vector.cpp
     */
    counting_vector & operator+=(bit_vector const & bit_vector)
    {
        for_each_set_bin(bit_vector,
                         [this](size_t const bin)
                         {
                             ++(*this)[bin];
                         });
        return *this;
    }

    /*!\brief Bin-wise subtracts the bits of a seqan::hibf::bit_vector.
     * \param bit_vector The seqan::hibf::bit_vector.
     * \attention The counting_vector must be at least as big as `bit_vector`.
     */
    counting_vector & operator-=(bit_vector const & bit_vector)
    {
        for_each_set_bin(bit_vector,
                         [this](size_t const bin)
                         {
                             assert((*this)[bin] > 0);
                             --(*this)[bin];
                         });
        return *this;
    }

    /*!\brief Bin-wise addition of two `seqan::hibf::counting_vector`s.
     * \param rhs The other seqan::hibf::counting_vector.
     * \attention The seqan::hibf::counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/counting_vector.cpp
     */
    counting_vector & operator+=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<value_t>());

        return *this;
    }

    /*!\brief Bin-wise substraction of two `seqan::hibf::counting_vector`s.
     * \param rhs The other seqan::hibf::counting_vector.
     * \attention The seqan::hibf::counting_vector must be at least as big as `rhs`.
     */
    counting_vector & operator-=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(),
                       this->end(),
                       rhs.begin(),
                       this->begin(),
                       [](auto a, auto b)
                       {
                           assert(a >= b);
                           return a - b;
                       });

        return *this;
    }

private:
    //!\brief Enumerates all bins of a seqan::hibf::bit_vector.
    template <typename on_bin_fn_t>
    void for_each_set_bin(bit_vector const & bit_vector, on_bin_fn_t && on_bin_fn)
    {
        assert(this->size() >= bit_vector.size()); // The counting vector may be bigger than what we need.
        size_t const words = (bit_vector.size() + 63u) >> 6;
        uint64_t const * const bitvector_raw = bit_vector.data();

        // Jump to the next 1 and return the number of places jumped in the bit_sequence
        auto jump_to_next_1bit = [](uint64_t & x)
        {
            auto const zeros = std::countr_zero(x);
            x >>= zeros; // skip number of zeros
            return zeros;
        };

        // Each iteration can handle 64 bits
        for (size_t batch = 0; batch < words; ++batch)
        {
            // get 64 bits starting at position `bit_pos`
            uint64_t bit_sequence = bitvector_raw[batch];

            // process each relative bin inside the bit_sequence
            for (size_t bin = batch << 6; bit_sequence != 0u; ++bin, bit_sequence >>= 1)
            {
                // Jump to the next 1 and
                bin += jump_to_next_1bit(bit_sequence);

                on_bin_fn(bin);
            }
        }
    }
};

} // namespace seqan::hibf
