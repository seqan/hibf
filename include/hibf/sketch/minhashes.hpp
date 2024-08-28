// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan::hibf::sketch::minhashes.
 */

#pragma once

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t, uint32_t
#include <span>    // for span
#include <vector>  // for vector

#include <cereal/access.hpp> // for access
#include <cereal/cereal.hpp> // for make_nvp, CEREAL_NVP

#include <hibf/misc/md_vector.hpp>
#include <hibf/platform.hpp>

namespace seqan::hibf::sketch
{

/*!\brief MinHash sketches design to be used for Locality sensitive hashing
 *
 * The partitioned HIBF as well as the Fast Layout need minHash sketches for a locality sentive hashing
 * algorithm that improves the user bin distribution.
 *
 * The sketches are specifically designed for this purpose:
 *
 * The `minhashes` struct keeps a table of several original MinHash sketches to be used for repeated
 * iterations of LSH.
 */
struct minhashes
{
    static constexpr uint64_t register_id_mask{15}; /// ...00001111
    static constexpr size_t num_sketches{16};
    static constexpr size_t sketch_size{40};

    //!\brief A table of sketches. For LSH we need multiple sketches, stored in a table.
    md_vector<uint64_t> table{}; // Each element (vector<uint64_t>) is a minhash.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    minhashes() = default;                              //!< Defaulted.
    minhashes(minhashes const &) = default;             //!< Defaulted.
    minhashes & operator=(minhashes const &) = default; //!< Defaulted.
    minhashes(minhashes &&) = default;                  //!< Defaulted.
    minhashes & operator=(minhashes &&) = default;      //!< Defaulted.
    ~minhashes() = default;                             //!< Defaulted.
    //!\brief construct from a vector of the smallest values in a set (sorted ascending).
    minhashes(std::vector<uint64_t> const & smallest_values);
    //!\}

    //!\brief Checks whether the minHash table is completely filled.
    bool is_valid() const;

    //!\brief Adds more minhash values to an existing but incomplete table.
    void fill_incomplete_sketches(std::span<uint64_t> const & more_smallest_values);

    //!\brief Pushes `value` to the heap if it is smaller than the current largest element.
    static void push_to_heap_if_smaller(uint64_t const value, std::vector<uint64_t> & heap);

private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t version{1};
        archive(CEREAL_NVP(version));

        // other members are const static currently
        archive(CEREAL_NVP(table));
    }
};

} // namespace seqan::hibf::sketch
