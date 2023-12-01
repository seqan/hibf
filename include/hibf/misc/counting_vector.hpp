// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::hibf::counting_vector.
 */

#pragma once

#include <algorithm>   // for transform
#include <bit>         // for countr_zero
#include <cassert>     // for assert
#include <cinttypes>   // for uint64_t
#include <climits>     // for CHAR_BIT
#include <concepts>    // for integral
#include <cstring>     // for size_t
#include <functional>  // for minus, plus
#include <type_traits> // for conditional_t
#include <vector>      // for vector

#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/misc/bit_vector.hpp>           // for bit_vector
#include <hibf/misc/divide_and_ceil.hpp>      // for divide_and_ceil
#include <hibf/misc/next_multiple_of_64.hpp>  // for next_multiple_of_64
#include <hibf/platform.hpp>                  // for HIBF_HAS_AVX512

#if HIBF_HAS_AVX512
#    include <simde/x86/avx512/add.h>   // for simde_mm512_add_epi16, simde_mm512_add_epi32, simde_mm512_add_...
#    include <simde/x86/avx512/load.h>  // for simde_mm512_load_si512
#    include <simde/x86/avx512/mov.h>   // for simde_mm512_maskz_mov_epi16, simde_mm512_maskz_mov_epi32, simd...
#    include <simde/x86/avx512/set1.h>  // for simde_mm512_set1_epi16, simde_mm512_set1_epi32, simde_mm512_se...
#    include <simde/x86/avx512/store.h> // for simde_mm512_store_si512
#    include <simde/x86/avx512/sub.h>   // for simde_mm512_sub_epi16, simde_mm512_sub_epi32, simde_mm512_sub_...
#    include <simde/x86/avx512/types.h> // for simde__m512i
#endif

namespace seqan::hibf
{

#if HIBF_HAS_AVX512
//!\cond
// Since the counting_vector can have different value types, we need specific SIMD instructions for each value type.
template <std::integral integral_t>
struct simd_mapping
{};

// CRTP base class for the simd_mapping, containg common functionality.
template <typename derived_t, std::integral integral_t>
    requires std::same_as<derived_t, simd_mapping<integral_t>>
struct simd_mapping_crtp
{
    // Let `B = sizeof(integral_t) * CHAR_BIT`, e.g. 8 for (u)int8_t, and 16 for (u)int16_t.
    // We can process `512 / B` bits of the bit_vector at once.
    static inline constexpr size_t bits_per_iterations = 512u / (sizeof(integral_t) * CHAR_BIT);
    // clang-format off
    // The type that is need to represent `bits_per_iterations` bits.
    using bits_type = std::conditional_t<bits_per_iterations == 64, uint64_t,
                      std::conditional_t<bits_per_iterations == 32, uint32_t,
                      std::conditional_t<bits_per_iterations == 16, uint16_t,
                      std::conditional_t<bits_per_iterations == 8, uint8_t, void>>>>;
    // clang-format on
    static_assert(!std::same_as<bits_type, void>, "Unsupported bits_type.");

    // Takes B bits from the bit_vector and expands them to a bits_type.
    // E.g., B = 8 : [1,0,1,1,0,0,1,0] -> [0...01, 0...00, 0...01, 0...01, 0...00, 0...00, 0...01, 0...00], where
    // each element is 64 bits wide.
    static inline constexpr auto expand_bits(bits_type const * const bits)
    {
        return derived_t::mm512_maskz_mov_epi(*bits, derived_t::mm512_set1_epi(1));
    }
};

// SIMD instructions for int8_t and uint8_t.
template <std::integral integral_t>
    requires (sizeof(integral_t) == 1)
struct simd_mapping<integral_t> : simd_mapping_crtp<simd_mapping<integral_t>, integral_t>
{
    static inline constexpr auto mm512_maskz_mov_epi = simde_mm512_maskz_mov_epi8;
    static inline constexpr auto mm512_set1_epi = simde_mm512_set1_epi8;
    static inline constexpr auto mm512_add_epi = simde_mm512_add_epi8;
    static inline constexpr auto mm512_sub_epi = simde_mm512_sub_epi8;
};

// SIMD instructions for int16_t and uint16_t.
template <std::integral integral_t>
    requires (sizeof(integral_t) == 2)
struct simd_mapping<integral_t> : simd_mapping_crtp<simd_mapping<integral_t>, integral_t>
{
    static inline constexpr auto mm512_maskz_mov_epi = simde_mm512_maskz_mov_epi16;
    static inline constexpr auto mm512_set1_epi = simde_mm512_set1_epi16;
    static inline constexpr auto mm512_add_epi = simde_mm512_add_epi16;
    static inline constexpr auto mm512_sub_epi = simde_mm512_sub_epi16;
};

// SIMD instructions for int32_t and uint32_t.
template <std::integral integral_t>
    requires (sizeof(integral_t) == 4)
struct simd_mapping<integral_t> : simd_mapping_crtp<simd_mapping<integral_t>, integral_t>
{
    static inline constexpr auto mm512_maskz_mov_epi = simde_mm512_maskz_mov_epi32;
    static inline constexpr auto mm512_set1_epi = simde_mm512_set1_epi32;
    static inline constexpr auto mm512_add_epi = simde_mm512_add_epi32;
    static inline constexpr auto mm512_sub_epi = simde_mm512_sub_epi32;
};

// SIMD instructions for int64_t and uint64_t.
template <std::integral integral_t>
    requires (sizeof(integral_t) == 8)
struct simd_mapping<integral_t> : simd_mapping_crtp<simd_mapping<integral_t>, integral_t>
{
    static inline constexpr auto mm512_maskz_mov_epi = simde_mm512_maskz_mov_epi64;
    static inline constexpr auto mm512_set1_epi = simde_mm512_set1_epi64;
    static inline constexpr auto mm512_add_epi = simde_mm512_add_epi64;
    static inline constexpr auto mm512_sub_epi = simde_mm512_sub_epi64;
};
//!\endcond
#endif

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
class counting_vector : public std::vector<value_t, seqan::hibf::contrib::aligned_allocator<value_t, 64u>>
{
private:
    //!\brief The base type.
    using base_t = std::vector<value_t, seqan::hibf::contrib::aligned_allocator<value_t, 64u>>;

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
     * \copydetails operator-=(bit_vector const &)
     * \details
     * ### Example
     *
     * \include test/snippet/ibf/counting_vector.cpp
     */
    counting_vector & operator+=(bit_vector const & bit_vector)
    {
        impl<operation::add>(bit_vector);
        return *this;
    }

    /*!\brief Bin-wise subtracts the bits of a seqan::hibf::bit_vector.
     * \param bit_vector The seqan::hibf::bit_vector.
     * \pre `counting_vector.size() >= bit_vector.size()` and either:
     *  * `bit_vector.size() % 64 == 0`
     * \pre or
     *  * `counting_vector.capacity() >= seqan::hibf::next_multiple_of_64(bit_vector.size())` and
     *    * &forall; bit &isin; [`bit_vector.size()`, `bit_vector.capacity()`) : `bit_vector[bit] == 0`<br>
     *      In practive, this condition might be violated by setting bits in a bit_vector and then shrinking it via
     *      \ref bit_vector::resize() "resize()".
     */
    counting_vector & operator-=(bit_vector const & bit_vector)
    {
        impl<operation::sub>(bit_vector);
        return *this;
    }

    /*!\brief Bin-wise addition of two `seqan::hibf::counting_vector`s.
     * \copydetails operator-=(counting_vector const &)
     * \details
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

    /*!\brief Bin-wise subtraction of two `seqan::hibf::counting_vector`s.
     * \param rhs The other seqan::hibf::counting_vector.
     * \pre `counting_vector.size() >= rhs.size()`
     */
    counting_vector & operator-=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::minus<value_t>());

        return *this;
    }

private:
    //!\brief The operation to perform.
    enum class operation
    {
        add,
        sub
    };

    //!\brief Bin-wise adds or subtracts the bits of a seqan::hibf::bit_vector.
    template <operation op>
    inline void impl(bit_vector const & bit_vector)
    {
        assert(this->size() >= bit_vector.size()); // The counting vector may be bigger than what we need.
#if HIBF_HAS_AVX512
        // AVX512BW: mm512_maskz_mov_epi, mm512_add_epi
        // AVX512F: mm512_set1_epi, _mm512_load_si512, _mm512_store_si512
        using simd = simd_mapping<value_t>;
        using bits_type = typename simd::bits_type;

        bits_type const * bit_vector_ptr = reinterpret_cast<bits_type const *>(bit_vector.data());
        value_t * counting_vector_ptr = base_t::data();

        size_t const bits = next_multiple_of_64(bit_vector.size());
        assert(bits <= this->capacity()); // Not enough memory reserved for AVX512 chunk access.
        size_t const iterations = bits / simd::bits_per_iterations;

        for (size_t iteration = 0; iteration < iterations;
             ++iteration, ++bit_vector_ptr, counting_vector_ptr += simd::bits_per_iterations)
        {
            simde__m512i load = simde_mm512_load_si512(counting_vector_ptr);
            if constexpr (op == operation::add)
            {
                load = simd::mm512_add_epi(load, simd::expand_bits(bit_vector_ptr));
            }
            else
            {
                load = simd::mm512_sub_epi(load, simd::expand_bits(bit_vector_ptr));
            }
            simde_mm512_store_si512(counting_vector_ptr, load);
        }
#else
        size_t const words = divide_and_ceil(bit_vector.size(), 64u);
        uint64_t const * const bit_vector_ptr = bit_vector.data();

        // Jump to the next 1 and return the number of jumped bits in value.
        auto jump_to_next_1bit = [](uint64_t & value)
        {
            auto const zeros = std::countr_zero(value);
            value >>= zeros; // skip number of zeros
            return zeros;
        };

        // Each iteration can handle 64 bits, i.e., one word.
        for (size_t iteration = 0; iteration < words; ++iteration)
        {
            uint64_t current_word = bit_vector_ptr[iteration];

            // For each set bit in the current word, add/subtract 1 to the corresponding bin.
            for (size_t bin = iteration * 64u; current_word != 0u; ++bin, current_word >>= 1)
            {
                // Jump to the next 1
                bin += jump_to_next_1bit(current_word);

                if constexpr (op == operation::add)
                {
                    ++(*this)[bin];
                }
                else
                {
                    --(*this)[bin];
                }
            }
        }
#endif
    }
};

} // namespace seqan::hibf
