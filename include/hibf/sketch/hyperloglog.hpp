// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-FileCopyrightText: 2013 Hideaki Ohno <hide.o.j55{at}gmail.com>
// SPDX-License-Identifier: BSD-3-Clause AND MIT

/*!\file
 * \author Felix Droop <felix.droop AT fu-berlin.de>
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::hibf::sketch::hyperloglog.
 */

#pragma once

#include <array>   // for array
#include <cstddef> // for size_t
#include <cstdint> // for uint64_t, uint8_t, uint32_t
#include <iosfwd>  // for istream, ostream
#include <vector>  // for vector

#include <cereal/access.hpp> // for access
#include <cereal/cereal.hpp> // for make_nvp, CEREAL_NVP

#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/platform.hpp>

namespace seqan::hibf::sketch
{

/*!\brief HyperLogLog estimates.
 * \ingroup hibf_sketch
 * \details
 * Original work by Hideaki Ohno. Major changes have been applied for bugfixes, 64-bit support, improvements, etc.
 * \see https://github.com/hideo55/cpp-HyperLogLog
 */
class hyperloglog
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    /*!\brief Default constructor.
     * \param[in] num_bits The bit width in [5,32].
     *
     * Allocates 2^`num_bits` bytes of memory.
     *
     * \throws std::invalid_argument if num_bits is not in [5,32].
     */
    hyperloglog(uint8_t const num_bits = 5u);
    hyperloglog(hyperloglog const &) = default;             //!< Defaulted.
    hyperloglog & operator=(hyperloglog const &) = default; //!< Defaulted.
    hyperloglog(hyperloglog &&) = default;                  //!< Defaulted.
    hyperloglog & operator=(hyperloglog &&) = default;      //!< Defaulted.
    ~hyperloglog() = default;                               //!< Defaulted.

    //!\}

    /*!\brief Adds a value.
     * \param[in] value The value to add.
     */
    void add(uint64_t const value);

    /*!\brief Estimates cardinality value.
     * \returns Estimated cardinality value.
     */
    double estimate() const;

    /*!\brief Merges another hyperloglog into this object.
     * \param[in] other The hyperloglog to be merged.
     * \details
     * This has the same effect as adding all values that were added to `other`.
     * \warning
     * Merging a hyperloglog with differing `bits` is undefined behaviour. In debug mode, this is an assertion instead.
     */
    void merge(hyperloglog const & other);

    /*!\brief Merges another hyperloglog and returns the new estimate.
     * \param[in] other The hyperloglog to be merged.
     * \returns Estimated cardinality value.
     * \details
     * \warning
     * Merging a hyperloglog with differing `bits` is undefined behaviour. In debug mode, this is an assertion instead.
     */
    double merge_and_estimate(hyperloglog const & other);

    /*!\brief Clears added values.
     * The size is unaffected.
     */
    void reset();

    /*!\brief Returns size of the internal data.
     * \returns Size in bytes.
     * The returned value is equivalent to 2^`bits`.
     */
    uint64_t data_size() const
    {
        return size;
    }

    /*!\brief Write the hyperloglog to a stream.
     * \param[in,out] os The output stream to write to.
     * \throws std::runtime_error if storing failed.
     */
    void store(std::ostream & os) const;

    /*!\brief Loads the hyperloglog from a stream.
     * \param[in] is The input stream where to read from.
     * \throws std::runtime_error if reading failed.
     */
    void load(std::istream & is);

private:
    //!\brief Used for estimation. Part of estimate E in the HyperLogLog publication.
    static constexpr std::array<float, 61> expectation_values = []() constexpr
    {
        std::array<float, 61> result{};
        for (size_t i = 0; i < 61; ++i)
            result[i] = 1.0f / (1ULL << i);
        return result;
    }();

    //!\brief The bit width. Also called precision, b, and p in other publications.
    uint8_t bits{};
    //!\brief Equivalent to 2^bits. Called m in original publication.
    uint64_t size{};
    //!\brief Mask used in add().
    uint64_t rank_mask{};
    //!\brief Equivalent to alpha * m^2.
    double normalization_factor{};
    //!\brief Internal data. Also called register in publications.
    std::vector<uint8_t, seqan::hibf::contrib::aligned_allocator<uint8_t, 32u>> data{};

    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t version{1};
        archive(CEREAL_NVP(version));

        archive(CEREAL_NVP(bits));
        archive(CEREAL_NVP(size));
        archive(CEREAL_NVP(rank_mask));
        archive(CEREAL_NVP(normalization_factor));
        archive(CEREAL_NVP(data));
    }
};

} // namespace seqan::hibf::sketch
