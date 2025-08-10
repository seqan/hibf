// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::hibf::interleaved_bloom_filter.
 */

#pragma once

#include <algorithm>   // for fill
#include <array>       // for array
#include <concepts>    // for integral, unsigned_integral, same_as
#include <cstdint>     // for uint16_t, uint64_t, uint32_t
#include <cstring>     // for size_t
#include <memory>      // for addressof
#include <ranges>      // for forward_range, range_value_t, range, range_reference_t, input...
#include <type_traits> // for remove_cvref_t, false_type
#include <vector>      // for vector, operator==

#include <cereal/cereal.hpp>           // for make_nvp
#include <cereal/macros.hpp>           // for CEREAL_SERIALIZE_FUNCTION_NAME
#include <cereal/types/base_class.hpp> // for base_class
#include <cereal/types/vector.hpp>     // IWYU pragma: keep

#include <hibf/cereal/concepts.hpp>            // for cereal_archive
#include <hibf/contrib/aligned_allocator.hpp>  // for aligned_allocator
#include <hibf/contrib/std/enumerate_view.hpp> // for enumerate
#include <hibf/misc/bit_vector.hpp>            // for bit_vector
#include <hibf/misc/counting_vector.hpp>       // for counting_vector
#include <hibf/misc/next_multiple_of_64.hpp>   // for next_multiple_of_64
#include <hibf/misc/unreachable.hpp>           // for unreachable, assert, forward
#include <hibf/platform.hpp>                   // for HIBF_CONSTEXPR_VECTOR, HIBF_HAS_AVX512

namespace seqan::hibf
{

// config.hpp -> misc/insert_iterator.hpp (Needs interleaved_bloom_filter to be a complete class)
struct config;

/*!\brief For testing: An accessor used to inspect private members of seqan::hibf::interleaved_bloom_filter.
 * \ingroup ibf
 */
struct inspector;

/*!\brief A strong type that represents the number of bins for the seqan::hibf::interleaved_bloom_filter.
 * \ingroup ibf
 * \qualifier strong
 */
struct bin_count
{
    size_t value;
};

/*!\brief A strong type that represents the number of bits for each bin in the seqan::hibf::interleaved_bloom_filter.
 * \ingroup ibf
 * \qualifier strong
 */
struct bin_size
{
    size_t value;
};

/*!\brief A strong type that represents the number of hash functions for the seqan::hibf::interleaved_bloom_filter.
 * \ingroup ibf
 * \qualifier strong
 */
struct hash_function_count
{
    size_t value;
};

/*!\brief A strong type that represents the bin index for the seqan::hibf::interleaved_bloom_filter.
 * \ingroup ibf
 * \qualifier strong
 */
struct bin_index
{
    size_t value;
};

/*!\brief The IBF binning directory. A data structure that efficiently answers set-membership queries for multiple bins.
 * \ingroup ibf
 * \implements seqan::hibf::cerealisable
 *
 * \details
 *
 * ### Binning Directory
 *
 * A binning directory is a data structure that can be used to determine set membership for elements.
 * For example, a common use case is dividing a database into a fixed number (e.g. 1024) bins by some means
 * of clustering (e.g. taxonomic binning or k-mer similarity clustering for genomic sequences).
 * For a query, the binning directory can now answer in which bins the query (probably) occurs.
 * In SeqAn we provide the Interleaved Bloom Filter (IBF) that can answer these queries efficiently.
 *
 * ### Interleaved Bloom Filter (IBF)
 *
 * The Interleaved Bloom Filter is a probabilistic data structure that extends the
 * [Bloom Filter](https://en.wikipedia.org/wiki/Bloom_filter).
 * A Bloom Filter can be thought of as a bitvector of length `n` and `h` hash functions and is used to determine set
 * membership. To insert data, the data is hashed by the `h` hash functions (returning values in `[0, n)`) and the
 * corresponding `h` positions in the bitvector are set to `1`. To query data, i.e. to determine whether the query
 * belongs to the set the Bloom Filter was built for, the query is hashed by the same `h` hash functions and the
 * corresponding positions are checked. If all `h` positions contain a `1`, the query is (probably) in the data set.
 * Since the Bloom Filter has variable length, the hashing is not bijective, i.e. it may return true for a set
 * membership query even though the query was never inserted into the Bloom Filter. Note that the Bloom Filter
 * will always return `true` if the query was inserted, i.e. there may be false positives, but no false negatives.
 *
 * The Interleaved Bloom Filter now applies the concept of a Bloom Filter to multiple sets and provides a *global*
 * data structure to determine set membership of a query in `b` data sets/bins.
 * Conceptually, a Bloom Filter is created for each bin using the same fixed length and fixed hash functions for each
 * filter. The resulting `b` Bloom Filters are then interleaved such that the `i`'th bit if each Bloom Filter are
 * adjacent to each other:
 * ```
 * Bloom Filter 0       Bloom Filter 1      Bloom Filter 2      Bloom Filter 3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * Where `x.y` denotes the `y`'th bit of the `x`'th Bloom Filter.
 * ```
 * Interleaved Bloom Filter
 * |0.0|1.0|2.0|3.0|0.1|1.1|2.1|3.1|0.2|1.2|2.2|3.2|0.3|1.3|2.3|3.3|
 * ```
 * A query can now be searched in all `b` bins by computing the `h` hash functions, retrieving the `h` sub-bitvectors of
 * length `b` starting at the positions indicated by the hash functions. The bitwise AND of these sub-bitvectors yields
 * the binningvector, a bitvector of length `b` where the `i`'th bit indicates set membership in the `i`'th bin.
 *
 * ### Querying
 *
 * | Agent               | Use case                                                                          | Agent's function                                                | Example Result   | Result type           |
 * |---------------------|-----------------------------------------------------------------------------------|-----------------------------------------------------------------|------------------|-----------------------|
 * | containment_agent() | Checks membership for a single value in all bins                                  | \ref containment_agent_type::bulk_contains() "bulk_contains()"  | `[0, 1, 1, 0]`   | bit_vector            |
 * | counting_agent()    | Counts occurrences for a range of values in all bins                              | \ref counting_agent_type::bulk_count() "bulk_count()"           | `[2, 23, 0, 7]`  | counting_vector       |
 * | membership_agent()  | Similar to counting, but returns all bin IDs that exceed (`>=`) a given threshold | \ref membership_agent_type::membership_for() "membership_for()" | `[1, 3]`         | std::vector<uint64_t> |
 *
 * ### Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * Additionally, concurrent calls to `emplace` are safe iff each thread handles a multiple of wordsize (=64) many bins.
 * For example, calls to `emplace` from multiple threads are safe if `thread_1` accesses bins 0-63, `thread_2` bins
 * 64-127, and so on.
 */
class interleaved_bloom_filter : private seqan::hibf::bit_vector
{
private:
    //!\brief The underlying datatype to use.
    using base_t = seqan::hibf::bit_vector;

    //!\brief Allow cereal to cast the IBF into its base class.
    template <typename t>
    friend struct cereal::base_class;

    //!\brief Allow access to private members. Used for util/inspect.
    friend struct seqan::hibf::inspector;

    //!\brief The number of bins specified by the user.
    size_t bins{};
    //!\brief The number of bins stored in the IBF (next multiple of 64 of `bins`).
    size_t technical_bins{};
    //!\brief The size of each bin in bits.
    size_t bin_size_{};
    //!\brief The number of bits to shift the hash value before doing multiplicative hashing.
    size_t hash_shift{};
    //!\brief The number of 64-bit integers needed to store `bins` many bits (e.g. `bins = 50` -> `bin_words = 1`).
    size_t bin_words{};
    //!\brief The number of hash functions.
    size_t hash_funs{};
    //!\brief Precalculated seeds for multiplicative hashing. We use large irrational numbers for a uniform hashing.
    static constexpr std::array<size_t, 5> hash_seeds{13572355802537770549ULL, // 2**64 / (e/2)
                                                      13043817825332782213ULL, // 2**64 / sqrt(2)
                                                      10650232656628343401ULL, // 2**64 / sqrt(3)
                                                      16499269484942379435ULL, // 2**64 / (sqrt(5)/2)
                                                      4893150838803335377ULL}; // 2**64 / (3*pi/5)

    /*!\brief Perturbs a value and fits it into the vector.
     * \param h The value to process.
     * \param seed The seed to use.
     * \returns A hashed value representing a position within the bounds of `data`.
     * \sa https://probablydance.com/2018/06/16/
     * \sa https://lemire.me/blog/2016/06/27
     */
    [[gnu::always_inline]] inline constexpr size_t hash_and_fit(size_t h, size_t const seed) const
    {
        h *= seed;
        assert(hash_shift < 64);
        h ^= h >> hash_shift;         // XOR and shift higher bits into lower bits
        h *= 11400714819323198485ULL; // = 2^64 / golden_ration, to expand h to 64 bit range
        // Use fastrange (integer modulo without division).
        h = static_cast<uint64_t>((static_cast<__uint128_t>(h) * static_cast<__uint128_t>(bin_size_)) >> 64);
        h *= technical_bins;
        return h;
    }

public:
    class containment_agent_type; // documented upon definition below
    template <std::integral value_t>
    class counting_agent_type;   // documented upon definition below
    class membership_agent_type; // documented upon definition below

    /*!\name Constructors, destructor and assignment
     * \{
     */
    interleaved_bloom_filter() = default;                                                 //!< Defaulted.
    interleaved_bloom_filter(interleaved_bloom_filter const &) = default;                 //!< Defaulted.
    interleaved_bloom_filter & operator=(interleaved_bloom_filter const &) = default;     //!< Defaulted.
    interleaved_bloom_filter(interleaved_bloom_filter &&) noexcept = default;             //!< Defaulted.
    interleaved_bloom_filter & operator=(interleaved_bloom_filter &&) noexcept = default; //!< Defaulted.
    ~interleaved_bloom_filter() = default;                                                //!< Defaulted.

    /*!\brief Construct an Interleaved Bloom Filter.
     * \param bins_ The number of bins.
     * \param size The bitvector size.
     * \param funs The number of hash functions. Default 2. At least 1, at most 5.
     * \param empty_bin_fraction The fraction of total technical bins that should be empty.
     *
     * \details
     *
     * Upon construction, `_bins` many bins are immediately accessible.
     * Additionally, enough space is reserve to accomodate a total of at least `_bins / empty_bin_fraction` many bins.
     * For example, with `bins_ = 64` and `empty_bin_fraction = 0.5`, the Interleaved Bloom Filter will reserve space
     * for 128 bins, 64 of which are immediately accessible and 64 of which are reserved for future use.
     *
     * ### Example
     *
     * \include test/snippet/ibf/interleaved_bloom_filter_constructor.cpp
     */
    interleaved_bloom_filter(seqan::hibf::bin_count const bins_,
                             seqan::hibf::bin_size const size,
                             seqan::hibf::hash_function_count const funs = seqan::hibf::hash_function_count{2u},
                             double const empty_bin_fraction = 0.0);

    /*!\brief Construct an Interleaved Bloom Filter.
     * \param configuration The seqan::hibf::config.
     * \param max_bin_elements Optional, the maximum number of unique elements in any bin.
     * \details
     *
     * If `max_bin_elements` is not passed, or `max_bin_elements` is 0, the maximum number of unique elements in any bin
     * will be determined automatically.
     *
     * `max_bin_elements` must be the maximum number of unique elements for any bin as evaluated with
     * `seqan::hibf::config::input_fn`.
     */
    interleaved_bloom_filter(config & configuration, size_t const max_bin_elements = 0u);
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Inserts a value into a specific bin.
     * \param[in] value The raw numeric value to process.
     * \param[in] bin The bin index to insert into.
     *
     * \details
     *
     * If `track_occupancy` is set to `true`, the occupancy of the bin is tracked.
     *
     * ### Example
     *
     * \include test/snippet/ibf/interleaved_bloom_filter_emplace.cpp
     */
    void emplace(size_t const value, bin_index const bin) noexcept;

    /*!\brief Clears a specific bin.
     * \param[in] bin The bin index to clear.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/interleaved_bloom_filter_clear.cpp
     */
    void clear(bin_index const bin) noexcept;

    /*!\brief Clears a range of bins.
     * \tparam rng_t The type of the range. Must model std::ranges::forward_range and the reference type must be
     *               seqan::hibf::bin_index.
     * \param[in] bin_range The range of bins to clear.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/interleaved_bloom_filter_clear.cpp
     */
    template <typename rng_t>
    void clear(rng_t && bin_range) noexcept
    {
        static_assert(std::ranges::forward_range<rng_t>, "The range of bins to clear must model a forward_range.");
        static_assert(std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>, bin_index>,
                      "The reference type of the range to clear must be seqan::hibf::bin_index.");
#ifndef NDEBUG
        for (auto && bin : bin_range)
            assert(bin.value < technical_bins);
#endif // NDEBUG

        for (size_t offset = 0, i = 0; i < bin_size_; offset += technical_bins, ++i)
            for (auto && bin : bin_range)
                (*this)[bin.value + offset] = 0;
    }

    /*!\brief Tries increasing the number of bins stored in the Interleaved Bloom Filter without reallocating memory.
     * \param[in] new_bin_count The new number of bins.
     * \returns `true` if the number of bins was set, `false` if the number of bins was not set.
     * \sa seqan::hibf::interleaved_bloom_filter::increase_bin_number_to
     *
     * \attention If the new bin count is greater than the old bin count and this function returns `true`, all
     * seqan::hibf::interleaved_bloom_filter::containment_agent_type and
     * seqan::hibf::interleaved_bloom_filter::counting_agent_type constructed for this Interleaved Bloom Filter are
     * invalidated.
     *
     * | Condition                         | Internal Condition     | Effect                        | Return value |
     * |-----------------------------------|------------------------|-------------------------------|--------------|
     * | New bin count > current bin count | resize is not required | bin count is set to new value | `true`       |
     * | New bin count > current bin count | resize is required     | none                          | `false`      |
     * | New bin count = current bin count | none                   | none                          | `true`       |
     * | New bin count < current bin count | none                   | none                          | `false`      |
     *
     * ### Example
     *
     * \include test/snippet/ibf/interleaved_bloom_filter_try_increase_bin_number_to.cpp
     */
    bool try_increase_bin_number_to(bin_count const new_bin_count) noexcept;

    /*!\brief Increases the number of bins stored in the Interleaved Bloom Filter.
     * \param[in] new_bin_count The new number of bins.
     * \throws std::invalid_argument If passed number of bins is smaller than current number of bins.
     *
     * \attention The new number of bins must be greater or equal to the current number of bins.
     * \attention This function invalidates all seqan::hibf::interleaved_bloom_filter::containment_agent_type and
     * seqan::hibf::interleaved_bloom_filter::counting_agent_type constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * The resulting `seqan::hibf::interleaved_bloom_filter` has an increased size proportional to the increase in the
     * `bin_words` (the number of 64-bit words needed to represent `bins` many bins), e.g.
     * resizing a `seqan::hibf::interleaved_bloom_filter` with 40 bins to 73 bins also increases the `bin_words` from 1 to
     * 2 and hence the new `seqan::hibf::interleaved_bloom_filter` will be twice the size.
     * This increase in size is necessary to avoid invalidating all computed hash functions.
     * If you want to add more bins while keeping the size constant, you need to rebuild the
     * `seqan::hibf::interleaved_bloom_filter`.
     *
     * ### Example
     *
     * \include test/snippet/ibf/interleaved_bloom_filter_increase_bin_number_to.cpp
     */
    void increase_bin_number_to(bin_count const new_bin_count);
    //!\}

    /*!\name Lookup
     * \{
     */
    /*!\brief Returns a seqan::hibf::interleaved_bloom_filter::containment_agent_type to be used for lookup.
     * \attention Calling seqan::hibf::interleaved_bloom_filter::increase_bin_number_to invalidates all
     * `seqan::hibf::interleaved_bloom_filter::containment_agent_type`s constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/containment_agent_construction.cpp
     * \sa seqan::hibf::interleaved_bloom_filter::containment_agent_type::bulk_contains
     */
    containment_agent_type containment_agent() const;

    /*!\brief Returns a seqan::hibf::interleaved_bloom_filter::counting_agent_type to be used for counting.
     * \attention Calling seqan::hibf::interleaved_bloom_filter::increase_bin_number_to invalidates all
     * `seqan::hibf::interleaved_bloom_filter::counting_agent_type`s constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/counting_agent_construction.cpp
     * \sa seqan::hibf::interleaved_bloom_filter::counting_agent_type::bulk_count
     */
    template <typename value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }

    /*!\brief Returns a seqan::hibf::interleaved_bloom_filter::membership_agent_type to be used for lookup.
     * \attention Calling seqan::hibf::interleaved_bloom_filter::increase_bin_number_to invalidates all
     * `seqan::hibf::interleaved_bloom_filter::membership_agent_type`s constructed for this Interleaved Bloom Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/containment_agent_construction.cpp
     * \sa seqan::hibf::interleaved_bloom_filter::membership_agent_type::membership_for
     */
    membership_agent_type membership_agent() const;
    //!\}

    /*!\name Capacity
     * \{
     */
    /*!\brief Returns the number of hash functions used in the Interleaved Bloom Filter.
     * \returns The number of hash functions.
     */
    size_t hash_function_count() const noexcept
    {
        return hash_funs;
    }

    /*!\brief Returns the number of bins that the Interleaved Bloom Filter manages.
     * \returns The number of bins.
     */
    size_t bin_count() const noexcept
    {
        return bins;
    }

    /*!\brief Returns the size of a single bin that the Interleaved Bloom Filter manages.
     * \returns The size in bits of a single bin.
     */
    size_t bin_size() const noexcept
    {
        return bin_size_;
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    size_t bit_size() const noexcept
    {
        return base_t::size();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    HIBF_CONSTEXPR_VECTOR bool operator==(interleaved_bloom_filter const &) const = default;
    //!\}

    /*!\name Access
     * \{
     */
    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A pointer to data of the underlying seqan::hibf::bit_vector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    using base_t::data;
    //!\}

    /*!\brief Contains the number of unique values inserted into each bin.
     * \details
     * Only contains non-zero values if `track_occupancy` is true.
     *
     * A value is unique if inserting it into the IBF would set at least one previously unset bit.
     */
    std::vector<size_t> occupancy{};

    //!\brief Whether to track the occupancy of the bins.
    bool track_occupancy{false};

    /*!\cond DEV
     * \brief The version of the HIBF.
     */
    static constexpr uint32_t version{1};

    /*!\brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan::hibf::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        uint32_t parsed_version{version};
        archive(cereal::make_nvp("version", parsed_version));
        archive(bins);
        archive(technical_bins);
        archive(bin_size_);
        archive(hash_shift);
        archive(bin_words);
        archive(hash_funs);
        archive(cereal::base_class<base_t>(this));
        archive(occupancy);
        archive(track_occupancy);
    }
    //!\endcond
};

/*!\brief Manages containment queries for the seqan::hibf::interleaved_bloom_filter.
 * \attention Calling seqan::hibf::interleaved_bloom_filter::increase_bin_number_to on `ibf` invalidates the
 * containment_agent.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/ibf/containment_agent_construction.cpp
 */
class interleaved_bloom_filter::containment_agent_type
{
private:
    //!\brief A pointer to the augmented seqan::hibf::interleaved_bloom_filter.
    interleaved_bloom_filter const * ibf_ptr{nullptr};

    //!\brief Stores access positions of augmented seqan::hibf::interleaved_bloom_filter.
    std::array<size_t, 5> bloom_filter_indices;

    //!\brief Stores the result of bulk_contains().
    bit_vector result_buffer;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    containment_agent_type() = default;                                               //!< Defaulted.
    containment_agent_type(containment_agent_type const &) = default;                 //!< Defaulted.
    containment_agent_type & operator=(containment_agent_type const &) = default;     //!< Defaulted.
    containment_agent_type(containment_agent_type &&) noexcept = default;             //!< Defaulted.
    containment_agent_type & operator=(containment_agent_type &&) noexcept = default; //!< Defaulted.
    ~containment_agent_type() = default;                                              //!< Defaulted.

    /*!\brief Construct a containment_agent_type from a seqan::hibf::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan::hibf::interleaved_bloom_filter.
     */
    explicit containment_agent_type(interleaved_bloom_filter const & ibf) :
        ibf_ptr(std::addressof(ibf)),
        result_buffer(ibf.bin_count())
    {}
    //!\}

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines containment (set membership) of a given value.
     * \param[in] value The raw value to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/containment_agent_bulk_contains.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan::hibf::interleaved_bloom_filter::containment_agent_type for each thread.
     */
    [[nodiscard]] bit_vector const & bulk_contains(size_t const value) & noexcept;

    //!\brief Deleted. Calling this function on a temporary would result in a dangling reference.
    [[nodiscard]] bit_vector const & bulk_contains(size_t const value) && noexcept = delete;
    //!\}
};

inline interleaved_bloom_filter::containment_agent_type interleaved_bloom_filter::containment_agent() const
{
    return interleaved_bloom_filter::containment_agent_type{*this};
}

/*!\brief Manages counting ranges of values for the seqan::hibf::interleaved_bloom_filter.
 * \attention Calling seqan::hibf::interleaved_bloom_filter::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/ibf/counting_agent.cpp
 */
template <std::integral value_t>
class interleaved_bloom_filter::counting_agent_type
{
private:
    //!\brief Store a seqan::hibf::interleaved_bloom_filter::containment_agent to call `bulk_contains`.
    containment_agent_type containment_agent;

    //!\brief Stores the result of bulk_count().
    counting_vector<value_t> result_buffer;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default;                                            //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default;                 //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default;     //!< Defaulted.
    counting_agent_type(counting_agent_type &&) noexcept = default;             //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) noexcept = default; //!< Defaulted.
    ~counting_agent_type() = default;                                           //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan::hibf::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan::hibf::interleaved_bloom_filter.
     */
    explicit counting_agent_type(interleaved_bloom_filter const & ibf) :
        containment_agent(ibf),
#if !HIBF_HAS_AVX512
        result_buffer(ibf.bin_count())
    {}
#else
        // AVX512 will access result_buffer's memory in chunks, so we need to make sure that we allocate enough memory
        // such that the last chunk is not out of bounds.
        result_buffer(next_multiple_of_64(ibf.bin_count())) // Ensure large enough capacity.
    {
        result_buffer.resize(ibf.bin_count()); // Resize to actual requested size.
                                               // Silences llvm's ASAN container-overflow warning.
#    if defined(_LIBCPP_VERSION)
#        if __has_feature(address_sanitizer)
#            if (_LIBCPP_VERSION < 200000)
        __sanitizer_annotate_contiguous_container(result_buffer.data(),
                                                  result_buffer.data() + result_buffer.capacity(),
                                                  result_buffer.data() + result_buffer.size(),
                                                  result_buffer.data() + result_buffer.capacity());
#            else
        std::__annotate_contiguous_container<allocator_type>(result_buffer.data(),
                                                             result_buffer.data() + result_buffer.capacity(),
                                                             result_buffer.data() + result_buffer.size(),
                                                             result_buffer.data() + result_buffer.capacity());
#            endif
#        endif
#    endif
    }
#endif

    //!\}

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::input_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan::hibf::interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector<value_t> const & bulk_count(value_range_t && values) & noexcept
    {
        static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, 0);

        for (auto && value : values)
            result_buffer += containment_agent.bulk_contains(value);

        return result_buffer;
    }

    //!\brief Deleted. Calling this function on a temporary would result in a dangling reference.
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector<value_t> const & bulk_count(value_range_t && values) && noexcept = delete;
    //!\}
};

/*!\brief Manages membership queries for the seqan::hibf::interleaved_bloom_filter.
 * \attention Calling seqan::hibf::interleaved_bloom_filter::increase_bin_number_to on `ibf` invalidates the
 * membership_agent.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/ibf/membership_agent_construction.cpp
 */
class interleaved_bloom_filter::membership_agent_type
{
private:
    //!\brief Store a seqan::hibf::interleaved_bloom_filter::counting_agent to call `bulk_count`.
    counting_agent_type<uint16_t> counting_agent;

    //!\brief Stores the result of membership_for().
    std::vector<uint64_t> result_buffer;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent_type() = default;                                              //!< Defaulted.
    membership_agent_type(membership_agent_type const &) = default;                 //!< Defaulted.
    membership_agent_type & operator=(membership_agent_type const &) = default;     //!< Defaulted.
    membership_agent_type(membership_agent_type &&) noexcept = default;             //!< Defaulted.
    membership_agent_type & operator=(membership_agent_type &&) noexcept = default; //!< Defaulted.
    ~membership_agent_type() = default;                                             //!< Defaulted.

    /*!\brief Construct a membership_agent_type from a seqan::hibf::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan::hibf::interleaved_bloom_filter.
     */
    explicit membership_agent_type(interleaved_bloom_filter const & ibf) :
        counting_agent(ibf),
        result_buffer(ibf.bin_count())
    {}
    //!\}

    /*!\brief Has no effect.
     * \details
     * This function exists to provide the same API as
     * seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent_type
     * \sa seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent_type::sort_results
     */
    void sort_results()
    {}

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership for all bins contained in this index, based on `values` and the `threshold`.
     * \param[in] values The values to process; must model std::ranges::forward_range.
     * \param[in] threshold Report a bin if there are at least this many hits.
     * \returns A vector of bin ids (index values) with successfull set membership query.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/ibf/membership_agent_membership_for.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan::hibf::interleaved_bloom_filter::membership_agent_type for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<uint64_t> const & membership_for(value_range_t && values,
                                                               uint16_t const threshold) & noexcept
    {
        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        auto & counts_per_bin = counting_agent.bulk_count(std::forward<value_range_t>(values));

        if (result_buffer.capacity() < counts_per_bin.size())
            seqan::hibf::unreachable();

        for (auto [bin, count] : seqan::stl::views::enumerate(counts_per_bin))
            if (count >= threshold)
                result_buffer.emplace_back(bin);

        return result_buffer;
    }

    //!\brief Deleted. Calling this function on a temporary would result in a dangling reference.
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<uint64_t> const & membership_for(value_range_t && values,
                                                               uint16_t const threshold) && noexcept = delete;
    //!\}

    //!\cond
    // GCOVR_EXCL_START
    template <typename t>
    struct always_false : std::false_type
    {};

    template <typename t = void>
    [[nodiscard]] bit_vector const & bulk_contains(size_t const) & noexcept
    {
        static_assert(
            always_false<t>::value,
            "The functionality of the membership_agent has been changed to be analogous to the "
            "hierarchical_interleaved_bloom_filter's version. As such, the interleaved_bloom_filter's membership_agent "
            "now only offers membership_for, which returns bin indices that exceed a given threshold when querying a "
            "range of values. bulk_contains is now part of a new agent: containment_agent. If you want to use "
            "bulk_contains, change ibf.membership_agent() to ibf.containment_agent().");
        seqan::hibf::unreachable();
    }

    [[nodiscard]] bit_vector const & bulk_contains(size_t const) && noexcept = delete;
    // GCOVR_EXCL_STOP
    //!\endcond
};

inline interleaved_bloom_filter::membership_agent_type interleaved_bloom_filter::membership_agent() const
{
    return interleaved_bloom_filter::membership_agent_type{*this};
}

} // namespace seqan::hibf
