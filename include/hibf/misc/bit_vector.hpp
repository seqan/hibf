// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::bit_vector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

// Modified by Enrico Seiler <enrico.seiler AT fu-berlin.de>:
// * To be a single header and non-template (was CRTP)
// * Changed default allocator to 64-byte aligned allocator
// * Changed `bit_vector operator~()` to be auto-vectorizable
// * Changed `binary_transform_impl` to be auto-vectorizable
// * Replaced `seqan3::detail::bits_of<chunk_type>` with `sizeof(chunk_type) * CHAR_BIT`
// * Changed `difference_type` to `size_type` in `operator[]`
// * Changed `1` to `1ULL` in `(1ULL << to_local_chunk_position(size()))` of `resize()`
// * Changed `bit_reference & operator=` to be in accordance with STL
// * Changed `resize` to do nothing when reducing size besides setting new size.

#pragma once

#include <algorithm>        // for __fn, for_each, all_of, any_of, copy, fill
#include <bit>              // for countr_zero
#include <cassert>          // for assert
#include <cinttypes>        // for uint64_t
#include <climits>          // for CHAR_BIT
#include <compare>          // for strong_ordering, operator==
#include <concepts>         // for assignable_from
#include <cstddef>          // for size_t, ptrdiff_t
#include <initializer_list> // for initializer_list
#include <iterator>         // for iter_reference_t, __fn, back_inserter, distance, iter_differen...
#include <memory>           // for allocator, assume_aligned, allocator_traits, __compressed_pair
#include <ranges>           // for __fn, begin, end
#include <stdexcept>        // for out_of_range
#include <string>           // for char_traits, operator+, to_string, operator""s
#include <type_traits>      // for conditional_t
#include <utility>          // for swap
#include <vector>           // for vector

#include <cereal/cereal.hpp>     // for binary_data
#include <cereal/macros.hpp>     // for CEREAL_LOAD_FUNCTION_NAME, CEREAL_SAVE_FUNCTION_NAME
#include <cereal/specialize.hpp> // for specialization, specialize

#include <hibf/cereal/concepts.hpp>           // for cereal_archive, cereal_text_archive
#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/platform.hpp>                  // for HIBF_CONSTEXPR_VECTOR, _LIBCPP_HAS_NO_ASAN, _LIBCPP_VERSION

namespace seqan::hibf
{

/*!\brief An bit vector.
 * \ingroup ibf
 *
 * \details
 *
 * Implements a bit vector on the basis of a std::vector with `uint64_t` as value type. The bit vector can be
 * dynamically resized and provides additional interfaces to apply efficient bit-operations on it.
 *
 * The reference type is a special proxy that provides access to a single bit. Note that it is not a real reference
 * but can be converter to a bool or assigned from a bool.
 */
class bit_vector :
    public std::vector<uint64_t,
                       typename std::allocator_traits<
                           seqan::hibf::contrib::aligned_allocator<bool, 64u>>::template rebind_alloc<uint64_t>>
{
private:
    //!\brief The allocator type.
    using allocator_t = seqan::hibf::contrib::aligned_allocator<bool, 64u>;

    //!\brief The base type.
    using base_t = std::vector<uint64_t, typename std::allocator_traits<allocator_t>::template rebind_alloc<uint64_t>>;

    //!\brief The type of the underlying chunk of bits.
    using chunk_type = uint64_t;

    /*!\brief The bit proxy returned as reference.
     *
     * \tparam is_const A bool that indicates a const proxy if the value is `true`, or a non-const proxy otherwise.
     *
     * \details
     *
     * This proxy is returned as a proxy for the reference type since a single bit cannot be addressed directly.
     * This proxy allows all operations that can be used on a single bit and is also implicitly convertible to a bool.
     * It cannot be default constructed and can only be instantiated with a particular bit from the bit vector and its
     * associated classes.
     */
    template <bool is_const>
    class bit_reference
    {
    private:
        //!\brief Befriend the bit vector so it can instantiate this proxy with a particular position.
        friend class bit_vector;

        //!\brief The const or non-const chunk type to be represented.
        using maybe_const_chunk_type = std::conditional_t<is_const, chunk_type const, chunk_type>;

        maybe_const_chunk_type * _chunk{}; //!< The underlying chunk.
        chunk_type _chunk_mask{};          //!< The mask selecting the bit for this chunk.

        /*!\brief Constructs the refernce with represented bit position within the container.
         *
         * \param[in] chunk A pointer to the chunk that contains the represented bit.
         * \param[in] local_chunk_position The position of the bit within the chunk.
         */
        constexpr bit_reference(maybe_const_chunk_type * chunk, size_type const local_chunk_position) noexcept :
            _chunk{chunk},
            _chunk_mask{static_cast<chunk_type>(1) << to_local_chunk_position(local_chunk_position)}
        {}

    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        bit_reference() = delete; //!< Deleted.
        bit_reference(bit_reference const & other) = default;
        bit_reference(bit_reference && other) = default;
        bit_reference & operator=(bit_reference const & other)
        {
            return *this = bool(other);
        }
        bit_reference & operator=(bit_reference && other) noexcept
        {
            return *this = bool(other);
        }

        /*!\brief Assigns a bit to the referenced bit.
         *
         * \param[in] bit The bit to set.
         */
        constexpr bit_reference & operator=(bool const bit) noexcept
        {
            bit ? set() : clear();
            return *this;
        }

        //!\overload
        // Needed to model std::output_iterator<bit_iterator, bool>, which requires the assignment to an const && version
        // of the proxy.
        constexpr bit_reference const & operator=(bool const bit) const noexcept
            requires (!is_const)
        {
            bit ? set() : clear();
            return *this;
        }
        //!\}

        //!\brief Converts this proxy to a bool.
        constexpr operator bool() const noexcept
        {
            return *_chunk & _chunk_mask;
        }

        //!\brief Flips the referenced bit.
        constexpr bit_reference & flip() noexcept
        {
            (*this) ? clear() : set();
            return *this;
        }

    private:
        //!\brief Sets the bit at the specific position.
        constexpr void set() noexcept
        {
            *_chunk |= _chunk_mask;
        }

        //!\brief Clears the bit at the specific position.
        constexpr void clear() noexcept
        {
            *_chunk &= ~_chunk_mask;
        }
    };

    /*!\brief A random access iterator over the bit vector.
     * \tparam is_const A bool that indicates a const iterator if the value is `true`, or a non-const iterator otherwise.
     */
    template <bool is_const>
    class bit_iterator
    {
    private:
        //!\brief Befriend the bit_iterator types with different constness.
        template <bool>
        friend class bit_iterator;

        //!\brief The type of the chunk.
        using maybe_const_chunk_type = std::conditional_t<is_const, chunk_type const, chunk_type>;

        maybe_const_chunk_type * _chunk{}; //!< The underlying chunk.
        size_type _chunk_position{};       //!< The bit position within the chunk.

    public:
        /*!\name Associated types
         * \{
         */
        using value_type = bool;                                   //!< The value type.
        using reference = bit_reference<is_const>;                 //!< The proxy type used as reference.
        using pointer = void;                                      //!\< The pointer type is void.
        using difference_type = std::ptrdiff_t;                    //!< The difference type.
        using iterator_category = std::random_access_iterator_tag; //!< The iterator category.
        using iterator_concept = std::random_access_iterator_tag;  //!< The iterator concept.
        //!\}

        /*!\name Constructors, destructor and assignment
         * \{
         */
        bit_iterator() = default; //!< Default.

        /*!\brief Constructs the iterator set to the begin of the given chunk.
         *
         * \param[in] chunk A pointer to the chunk that contains the represented bit.
         */
        explicit constexpr bit_iterator(maybe_const_chunk_type * chunk) noexcept : _chunk{chunk}, _chunk_position{0}
        {}

        /*!\brief Copies from a non-const iterator.
         *
         * \param[in] other The other non-const iterator to copy from.
         */
        constexpr bit_iterator(bit_iterator<!is_const> const & other) noexcept
            requires (is_const)
            : _chunk{other._chunk}, _chunk_position{other._chunk_position}
        {}
        //!\}

        /*!\name Element access
         * \{
         */
        //!\brief Returns the currently pointer-to element.
        constexpr reference operator*() const noexcept
        {
            return reference{_chunk, _chunk_position};
        }

        //!\brief Returns the element `count` positions away from the current iterator position.
        constexpr reference operator[](difference_type const count) const noexcept
        {
            return *((*this) + count);
        }
        //!\}

        /*!\name Arithmetic operator
         * \{
         */
        //!\brief Increments the iterator by one.
        constexpr bit_iterator & operator++() noexcept
        {
            _chunk += !static_cast<bool>(to_local_chunk_position(++_chunk_position));
            return *this;
        }

        //!\brief Increments the iterator by one and returns the iterator before the increment.
        constexpr bit_iterator operator++(int) noexcept
        {
            bit_iterator tmp{*this};
            ++(*this);
            return tmp;
        }

        //!\brief Advances the iterator by `count` many elements.
        constexpr bit_iterator & operator+=(difference_type const count) noexcept
        {
            //           chunk:|    0   |    1   |    2   |    3   |    4   |    5   |
            //                 |--------|--------|--------|--------|-x------|--------|
            //  chunk_position:|01234567|01234567|01234567|01234567|01234567|01234567|
            // global position:|01234567|89012345|67890123|45678901|23456789|01234567|
            //                 |0       |  1     |    2   |      3 |        |4       |
            if (count < 0)
            {
                size_type updated_count = modulo_mask - to_local_chunk_position(_chunk_position) - count;
                _chunk_position = modulo_mask - to_local_chunk_position(updated_count);
                _chunk -= to_chunk_position(updated_count);
                //(to_chunk_position(-count) + (old_chunk_position < _chunk_position));
            }
            else
            {
                _chunk += to_chunk_position(to_local_chunk_position(_chunk_position) + count);
                _chunk_position = to_local_chunk_position(_chunk_position + count);
            }

            return *this;
        }

        //!\brief Returns a new iterator advanced by `count` many elements.
        constexpr bit_iterator operator+(difference_type const count) const noexcept
        {
            bit_iterator tmp{*this};
            return tmp += count;
        }

        //!\brief Returns a new iterator advanced by `count` many elements.
        friend constexpr bit_iterator operator+(difference_type const count, bit_iterator rhs) noexcept
        {
            return rhs + count;
        }

        //!\brief Decrements the iterator by one.
        constexpr bit_iterator & operator--() noexcept
        {
            _chunk -= !static_cast<bool>(to_local_chunk_position(--_chunk_position));
            return *this;
        }

        //!\brief Decrements the iterator by one and returns the iterator before the decrement.
        constexpr bit_iterator operator--(int) noexcept
        {
            bit_iterator tmp{*this};
            --(*this);
            return tmp;
        }

        //!\brief Advances the iterator by `count` many elements.
        constexpr bit_iterator & operator-=(difference_type const count) noexcept
        {
            return *this += -count;
        }

        //!\brief Returns a new iterator advances by `count` many elements.
        constexpr bit_iterator operator-(difference_type const count) const noexcept
        {
            bit_iterator tmp{*this};
            return tmp -= count;
        }

        //!\brief Returns the distance between `this` and the `rhs` iterator.
        template <bool is_const_other>
        constexpr difference_type operator-(bit_iterator<is_const_other> rhs) const noexcept
        {
            return ((_chunk - rhs._chunk) << division_mask) -     // number of bits between chunks.
                   to_local_chunk_position(rhs._chunk_position) + // minus the first bits in rhs.
                   to_local_chunk_position(_chunk_position);      // plus the first bits of the lhs
        }
        //!\}

        /*!\name Comparison operators
         * \{
         */
        //!\brief Compares with another iterator.
        template <bool is_const_other>
        bool operator==(bit_iterator<is_const_other> const & rhs) const
        {
            return _chunk == rhs._chunk
                && (to_local_chunk_position(_chunk_position) == to_local_chunk_position(rhs._chunk_position));
        }

        //!\brief Compares the two iterator by their chunk position and local chunk position.
        template <bool is_const_other>
        std::strong_ordering operator<=>(bit_iterator<is_const_other> const & rhs) const
        {
            if (std::strong_ordering order = _chunk <=> rhs._chunk; order == std::strong_ordering::equivalent)
                return to_local_chunk_position(_chunk_position) <=> to_local_chunk_position(rhs._chunk_position);
            else
                return order;
        }
        //!\}
    };

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The iterator over the bits.
    using iterator = bit_iterator<false>;
    //!\brief The const iterator over the bits.
    using const_iterator = bit_iterator<true>;
    //!\brief The value type is `bool`.
    using value_type = std::iter_value_t<iterator>;
    //!\brief The reference type which is implemented as a proxy.
    using reference = std::iter_reference_t<iterator>;
    //!\brief The const reference type which is implemented as a proxy.
    using const_reference = std::iter_reference_t<const_iterator>;
    //!\brief The size_type.
    using size_type = size_t;
    //!\brief The difference type.
    using difference_type = std::iter_difference_t<iterator>;
    //!\brief The allocator type to use.
    using allocator_type = allocator_t;
    //!\}

private:
    //!\brief The number of bits represented in one chunk, e.g. 64.
    static constexpr size_type chunk_size = sizeof(chunk_type) * CHAR_BIT;
    //!\brief The mask used for the modulo operations using the bitwise and operator, e.g. & 63.
    static constexpr size_type modulo_mask = chunk_size - 1u;
    //!\brief The mask used for the division operations using bitwise shift operator, e.g. >> 6.
    static constexpr size_type division_mask = std::countr_zero(chunk_size);

    //!\brief The number of elements.
    size_type _size{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief The default constructor which optionally sets the allocator.
    HIBF_CONSTEXPR_VECTOR bit_vector(allocator_type const & alloc = allocator_type{}) : base_t{alloc}
    {}
    HIBF_CONSTEXPR_VECTOR bit_vector(bit_vector const &) = default;             //!< Default.
    HIBF_CONSTEXPR_VECTOR bit_vector(bit_vector &&) = default;                  //!< Default.
    HIBF_CONSTEXPR_VECTOR bit_vector & operator=(bit_vector const &) = default; //!< Default.
    HIBF_CONSTEXPR_VECTOR bit_vector & operator=(bit_vector &&) = default;      //!< Default.
    HIBF_CONSTEXPR_VECTOR ~bit_vector() = default;                              //!< Default.

    /*!\brief Constructs the bit vector with `count` copies of elements with value `bit`.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] bit The bit to set during initialisation.
     * \param[in] alloc The allocator to use [optional].
     */
    HIBF_CONSTEXPR_VECTOR
    bit_vector(size_type const count, bool const bit, allocator_type const & alloc = allocator_type{}) : base_t{alloc}
    {
        assign(count, bit);
    }

    /*!\brief Constructs the container initialised with the elements in `list`.
     *
     * \param[in] list An initialiser list with the bits set.
     * \param[in] alloc The allocator to use [optional].
     */
    HIBF_CONSTEXPR_VECTOR bit_vector(std::initializer_list<bool> list,
                                     allocator_type const & alloc = allocator_type{}) :
        base_t{alloc}
    {
        assign(list);
    }

    /*!\brief Constructs the container with `count` default-inserted instances of `bool`. No copies are made.
     *
     * \param[in] count The number of elements to create the bit vector with.
     * \param[in] alloc The allocator to use [optional].
     */
    HIBF_CONSTEXPR_VECTOR bit_vector(size_type const count, allocator_type const & alloc = allocator_type{}) :
        bit_vector{count, bool{}, alloc}
    {}
    //!\}

    /*!\name Member functions
     * \{
     */

    /*!\brief Assigns values to the container.
     *
     * \tparam iterator_t The type of the iterator; must model std::input_iterator.
     * \tparam sentinel_t The type of the sentinel; must model std::sentinel_for `iterator_t`.
     *
     * \param[in] first An first element to copy elements from.
     * \param[in] last The end of the range to copy elements from.
     *
     * \details
     *
     * Replaces the contents with copies of the range `[first, last)'. The behaviour is undefined if either argument
     * is an iterator to `*this`.
     *
     * All iterators, pointers and references to the elements of the container are invalidated.
     * The past-the-end iterator is also invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in distance between first and last.
     */
    template <std::input_iterator iterator_t, std::sentinel_for<iterator_t> sentinel_t>
        requires std::assignable_from<bool &, std::iter_reference_t<iterator_t>>
    constexpr void assign(iterator_t first, sentinel_t last)
    {
        bit_vector tmp{}; // To ensure strong exception guarantee.
        if constexpr (std::sized_sentinel_for<sentinel_t, iterator_t>)
            tmp.reserve(std::ranges::distance(first, last));

        std::ranges::copy(first, last, std::back_inserter(tmp));

        // ----- no exception after this.
        swap(tmp);
        set_new_size(std::ranges::distance(begin(), end()));
    }

    /*!\brief Assigns values to the container.
     *
     * \param[in] ilist The initialiser list with the elements to insert.
     *
     * \details
     *
     * Replaces the contents with the elements from the initializer list ilist.
     *
     * All iterators, pointers and references to the elements of the container are invalidated.
     * The past-the-end iterator is also invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in `ilist.size()`.
     */
    constexpr void assign(std::initializer_list<bool> const & ilist)
    {
        assign(std::ranges::begin(ilist), std::ranges::end(ilist));
    }

    /*!\brief Assigns values to the container.
     *
     * \param[in] count The new size of the container.
     * \param[in] bit The value to initialize elements of the container with.
     *
     * \details
     *
     * Replaces the contents with `count` copies of value `bit`.
     *
     * All iterators, pointers and references to the elements of the container are invalidated.
     * The past-the-end iterator is also invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in count.
     */
    HIBF_CONSTEXPR_VECTOR void assign(size_type const count, bool const bit)
    {
        resize(count, bit);
        std::ranges::for_each(*as_base(),
                              [value = fill_chunk(bit)](chunk_type & chunk)
                              {
                                  chunk = value;
                              });
    }
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Access specified element.
    HIBF_CONSTEXPR_VECTOR reference operator[](size_type const position) noexcept
    {
        assert(position < size());

        return *std::ranges::next(begin(), position);
    }

    //!\brief Access specified element.
    HIBF_CONSTEXPR_VECTOR const_reference operator[](size_type const position) const noexcept
    {
        assert(position < size());

        return *std::ranges::next(begin(), position);
    }

    /*!\brief Access the last element.
     *
     * \returns A reference to the last element in the container.
     *
     * \details
     *
     * Calling back on an empty container causes underfined behaviour.
     *
     * ### Exception
     *
     * Throws nothing.
     *
     * ### Complexity
     *
     * Constant.
     */
    HIBF_CONSTEXPR_VECTOR reference back() noexcept
    {
        assert(!empty()); // Calling on empty container is undefined behaviour.

        return (*this)[size() - 1u];
    }

    //!\overload
    HIBF_CONSTEXPR_VECTOR const_reference back() const noexcept
    {
        assert(!empty()); // Calling on empty container is undefined behaviour.

        return (*this)[size() - 1u];
    }

    //!\brief Checks if all bits are set to `true`.
    constexpr bool all() const noexcept
    {
        constexpr chunk_type mask = ~static_cast<chunk_type>(0);
        return std::ranges::all_of(*as_base(),
                                   [](chunk_type const & chunk)
                                   {
                                       return chunk == mask;
                                   });
    }

    //!\brief Checks if any bit is set to `true`.
    constexpr bool any() const noexcept
    {
        constexpr chunk_type mask = static_cast<chunk_type>(0);
        return std::ranges::any_of(*as_base(),
                                   [](chunk_type const & chunk)
                                   {
                                       return chunk | mask;
                                   });
    }

    //!\brief Checks if none of the bits is set to `true`.
    constexpr bool none() const noexcept
    {
        return !any();
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Returns the number of elements.
    constexpr size_type size() const noexcept
    {
        return _size;
    }

    //!\brief Checks wether the container is empty.
    constexpr bool empty() const noexcept
    {
        return _size == 0u;
    }

    //!\brief Returns the capacity.
    HIBF_CONSTEXPR_VECTOR size_type capacity() const noexcept
    {
        return base_t::capacity() * chunk_size;
    }

    /*!\brief Reserves storage.
     *
     * \param[in] new_capacity The new capacity of the bit vector.
     *
     * \details
     *
     * Increase the capacity of the vector to a value that's greater or equal to new_capacity. If new_capacity is
     * greater than the current capacity(), new storage is allocated, otherwise the method does nothing.
     * reserve() does not change the size of the vector. If new_capacity is greater than capacity(), all iterators,
     * including the past-the-end iterator, and all references to the elements are invalidated. Otherwise, no
     * iterators or references are invalidated.
     *
     * ### Exceptions
     *
     * std::length_error if `new_capacity > max_size()` or any exception thrown by allocator_t::allocate().
     * If an exception is thrown this function has no effect (strong exception guarantee).
     */
    HIBF_CONSTEXPR_VECTOR void reserve(size_type const new_capacity)
    {
        base_t::reserve(host_size_impl(new_capacity));
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Adds an element to the end.
     *
     * \param bit The bit to add to the end.
     *
     * \details
     *
     * Appends the given element value to the end of the container.
     *
     * If the new size() is greater than capacity() then all iterators and references
     * (including the past-the-end iterator) are invalidated.
     * Otherwise only the past-the-end iterator is invalidated.
     *
     * ### Exception
     *
     * If an exception is thrown (which can be due to allocator_t::allocate(), this function has no effect
     * (strong exception guarantee).
     *
     * ### Complexity
     *
     * Amortised constant.
     */
    HIBF_CONSTEXPR_VECTOR void push_back(bool bit)
    {
        size_t const new_size = size() + 1u;
        resize(new_size);
        // ---- no exception after this point.
        set_new_size(new_size);
        back() = bit; // set the bit.
    }

    //!\brief Changes the number of elements stored, where additional copies of `bit` are appended.
    HIBF_CONSTEXPR_VECTOR void resize(size_type const count, bool const bit = {})
    {
        base_t::resize(host_size_impl(count));

        size_t const old_size = size();
        set_new_size(count);

        // If bit is true and we increase the size.
        if (bit && size() > old_size)
            std::ranges::fill(begin() + old_size, end(), bit);
    }

    //!\brief Erases all elements. After this call, size() returns zero. capacity() remains unchanged.
    HIBF_CONSTEXPR_VECTOR void clear() noexcept
    {
        base_t::clear();
        _size = 0u;
    }

    //!\brief Performs binary AND between `this` and `rhs`.
    constexpr bit_vector & operator&=(bit_vector const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return binary_transform_impl(rhs,
                                     [](auto const & left_chunk, auto const & right_chunk)
                                     {
                                         return left_chunk & right_chunk;
                                     });
    }

    //!\brief Performs binary OR between `this` and `rhs`.
    constexpr bit_vector & operator|=(bit_vector const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return binary_transform_impl(rhs,
                                     [](auto const & left_chunk, auto const & right_chunk)
                                     {
                                         return left_chunk | right_chunk;
                                     });
    }

    //!\brief Performs binary XOR between `this` and `rhs`.
    constexpr bit_vector & operator^=(bit_vector const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return binary_transform_impl(rhs,
                                     [](auto const & left_chunk, auto const & right_chunk)
                                     {
                                         return left_chunk ^ right_chunk;
                                     });
    }

    //!\brief Performs binary NOT.
    HIBF_CONSTEXPR_VECTOR bit_vector operator~() const noexcept
    {
        bit_vector tmp(size());

        tmp.binary_transform_impl(*this,
                                  [](auto const &, auto const & right_chunk)
                                  {
                                      return ~right_chunk;
                                  });

        return tmp;
    }

    //!\brief Performs binary AND.
    HIBF_CONSTEXPR_VECTOR friend bit_vector operator&(bit_vector lhs, bit_vector const & rhs) noexcept
    {
        return lhs &= rhs;
    }

    //!\brief Performs binary OR.
    HIBF_CONSTEXPR_VECTOR friend bit_vector operator|(bit_vector lhs, bit_vector const & rhs) noexcept
    {
        return lhs |= rhs;
    }

    //!\brief Performs binary XOR.
    HIBF_CONSTEXPR_VECTOR friend bit_vector operator^(bit_vector lhs, bit_vector const & rhs) noexcept
    {
        return lhs ^= rhs;
    }

    //!\brief Computes the bitwise `a &= ~b` operator without an additional copy.
    constexpr bit_vector & and_not(bit_vector const & rhs) noexcept
    {
        assert(rhs.size() == size());

        return binary_transform_impl(rhs,
                                     [](auto const & left_chunk, auto const & right_chunk)
                                     {
                                         return left_chunk & ~right_chunk;
                                     });
    }

    //!\brief Flips all bits in-place.
    constexpr bit_vector & flip() noexcept
    {
        std::ranges::for_each(*as_base(),
                              [](chunk_type & chunk)
                              {
                                  chunk = ~chunk;
                              });
        return *this;
    }

    //!\brief Flips the bit at the given position.
    HIBF_CONSTEXPR_VECTOR bit_vector & flip(size_type position)
    {
        using namespace std::literals;

        if (position >= size())
            throw std::out_of_range{"The given posisiton "s + std::to_string(position) + " is out of the range [0, "s
                                    + std::to_string(size()) + ")!"s};

        (*this)[position].flip();
        return *this;
    }

    //!\brief Exchanges the contents of the container with those of others.
    HIBF_CONSTEXPR_VECTOR void swap(bit_vector & other) noexcept
    {
        base_t::swap(*other.as_base());
        std::swap(_size, other._size);
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator to the beginning.
    HIBF_CONSTEXPR_VECTOR iterator begin() noexcept
    {
        return iterator{base_t::data()};
    }

    //!\overload
    HIBF_CONSTEXPR_VECTOR const_iterator begin() const noexcept
    {
        return const_iterator{base_t::data()};
    }

    //!\overload
    HIBF_CONSTEXPR_VECTOR const_iterator cbegin() const noexcept
    {
        return begin();
    }

    //!\brief Returns an iterator to the end.
    HIBF_CONSTEXPR_VECTOR iterator end() noexcept
    {
        return begin() + size();
    }

    //!\overload
    HIBF_CONSTEXPR_VECTOR const_iterator end() const noexcept
    {
        return begin() + size();
    }

    //!\overload
    HIBF_CONSTEXPR_VECTOR const_iterator cend() const noexcept
    {
        return end();
    }
    //!\}

    [[gnu::always_inline]] inline HIBF_CONSTEXPR_VECTOR chunk_type * data() noexcept
    {
        return std::assume_aligned<allocator_type::alignment>(base_t::data());
    }

    [[gnu::always_inline]] inline HIBF_CONSTEXPR_VECTOR chunk_type const * data() const noexcept
    {
        return std::assume_aligned<allocator_type::alignment>(base_t::data());
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan::hibf::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly.
     */
    template <cereal_archive archive_t>
    void CEREAL_LOAD_FUNCTION_NAME(archive_t & archive)
    {
        // Not using `cereal::make_size_tag(_size)`, because the size tag is inferred for text (XML/JSON) archives.
        // For text archives, `cereal::make_size_tag(_size)` would be the number of elements serialised in the for-loop.
        // E.g., `_size == 100` would store `2` (`== host_size_impl(_size)`).
        archive(_size);
        size_t const vector_size = host_size_impl(_size);

        resize_for_overwrite(vector_size);

        if constexpr (cereal_text_archive<archive_t>)
        {
            for (auto && v : *as_base())
                archive(v);
        }
        else
        {
            archive(cereal::binary_data(data(), vector_size * sizeof(chunk_type)));
        }
    }

    //!\copydoc load
    template <cereal_archive archive_t>
    void CEREAL_SAVE_FUNCTION_NAME(archive_t & archive) const
    {
        // Not using `cereal::make_size_tag(_size)`, because the size tag is inferred for text (XML/JSON) archives.
        // For text archives, `cereal::make_size_tag(_size)` would be the number of elements serialised in the for-loop.
        // E.g., `_size == 100` would store `2` (`== host_size_impl(_size)`).
        archive(_size);

        if constexpr (cereal_text_archive<archive_t>)
        {
            for (auto && v : *as_base())
                archive(v);
        }
        else
        {
            archive(cereal::binary_data(data(), base_t::size() * sizeof(chunk_type)));
        }
    }
    //!\endcond

private:
// If nothing else works: Just use `resize`.
#ifndef HIBF_UNINITIALISED_RESIZE
    HIBF_CONSTEXPR_VECTOR inline void resize_for_overwrite(size_t const size)
    {
        base_t::resize(size);
    }
#else
// libc++: reinterpret_cast to a struct that has the same layout as std::vector.
// All internal members are private.
#    if defined(_LIBCPP_VERSION)
    inline void resize_for_overwrite(size_t const size)
    {
        struct fake_vector
        {
            using allocator_t = typename base_t::allocator_type;
            using pointer = typename std::allocator_traits<allocator_t>::pointer;

            pointer begin;
            pointer end;
            std::__compressed_pair<pointer, allocator_t> end_cap;
        };

        static_assert(sizeof(fake_vector) == sizeof(base_t));
        static_assert(alignof(fake_vector) == alignof(base_t));

        if (size > base_t::capacity())
            base_t::reserve(size);

// Annotate the new memory as contiguous container for llvm's address sanitizer.
#        ifndef _LIBCPP_HAS_NO_ASAN
        __sanitizer_annotate_contiguous_container(base_t::data(),
                                                  base_t::data() + base_t::capacity(),
                                                  base_t::data() + base_t::size(),
                                                  base_t::data() + size);
#        endif

        fake_vector & vec = reinterpret_cast<fake_vector &>(*this);
        vec.end = vec.begin + size;
    }
// libstdc++: The internal members are protected, so we can access them.
#    else
    HIBF_CONSTEXPR_VECTOR inline void resize_for_overwrite(size_t const size)
    {
        if (size > base_t::capacity())
            base_t::reserve(size);

        this->_M_impl._M_finish = this->_M_impl._M_start + size;
    }
#    endif
#endif

    //!\brief Performs the binary bitwise-operation on the underlying chunks.
    template <typename binary_operator_t>
    constexpr bit_vector & binary_transform_impl(bit_vector const & rhs, binary_operator_t && op) noexcept
    {
        chunk_type * const lhs_data = data();
        chunk_type const * const rhs_data = rhs.data();
        size_type const size = host_size_impl(this->size());

        for (size_t i = 0; i < size; ++i)
            lhs_data[i] = op(lhs_data[i], rhs_data[i]);

        return *this;
    }

    //!\brief Computes the minimal size needed for the host vector.
    //!\param[in] count The number of bits to allocate memory for.
    constexpr size_type host_size_impl(size_type const count) const noexcept
    {
        return chunks_needed(count);
    }

    //!\brief Sets the new size.
    constexpr void set_new_size(size_type const new_size) noexcept
    {
        _size = new_size;
    }

    //!\brief Casts `this` to its base class.
    constexpr base_t const * as_base() const noexcept
    {
        return static_cast<base_t const *>(this);
    }

    //!\overload
    constexpr base_t * as_base() noexcept
    {
        return static_cast<base_t *>(this);
    }

    //!\brief Returns how many chunks are needed to store `count` many elements.
    constexpr size_type chunks_needed(size_type const count) const noexcept
    {
        return (count + 63u) >> 6; // ceil(count/64)
    }

    //!\brief Returns a new chunk filled with the given bit.
    constexpr chunk_type fill_chunk(bool const bit) const noexcept
    {
        return (bit) ? ~chunk_type{} : chunk_type{};
    }

    //!\brief Converts the position to the local position within the chunk.
    static constexpr size_type to_local_chunk_position(size_type const position) noexcept
    {
        return position & modulo_mask; // e.g. position % 64
    }

    //!\brief Converts the position to the chunk position.
    static constexpr size_type to_chunk_position(size_type const position) noexcept
    {
        return position >> division_mask; // e.g. position / 64
    }
};

} // namespace seqan::hibf

//!\cond
// See https://uscilab.github.io/cereal/serialization_functions.html#inheritance
// seqan::hibf::bit_vector's base class is std::vector
// If we include <cereal/types/vector.hpp> for std::vector serialisation (e.g., HIBF),
// cereal provides these as non-member load/save functions.
// Since both load/save non-member functions (std::vector) and load/save member functions (bit_vector) are available,
// cereal needs to be told which one to use.
namespace cereal
{

template <typename archive_t>
struct specialize<archive_t, seqan::hibf::bit_vector, cereal::specialization::member_load_save>
{};

} // namespace cereal
//!\endcond
