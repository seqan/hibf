// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for ptrdiff_t
#include <functional> // for equal_to
#include <iterator>   // for output_iterator_tag
#include <memory>     // for addressof
#include <vector>     // for vector

#include <hibf/contrib/robin_hood.hpp> // for unordered_flat_set, hash
#include <hibf/platform.hpp>
#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

// IWYU pragma: private, include <hibf/config.hpp>

namespace seqan::hibf
{

class insert_iterator
{
public:
    using iterator_category = std::output_iterator_tag;
    using value_type = void;
    using difference_type = ptrdiff_t;
    using pointer = void;
    using reference = void;

    insert_iterator() = delete;
    insert_iterator(insert_iterator const &) = default;
    insert_iterator(insert_iterator &&) = default;
    insert_iterator & operator=(insert_iterator const &) = default;
    insert_iterator & operator=(insert_iterator &&) = default;
    ~insert_iterator() = default;

    explicit constexpr insert_iterator(robin_hood::unordered_flat_set<uint64_t> & set) :
        set{std::addressof(set)},
        type{data_type::unordered_set}
    {}

    explicit constexpr insert_iterator(std::vector<uint64_t> & vec) : vec{std::addressof(vec)}, type{data_type::vector}
    {}

    explicit constexpr insert_iterator(sketch::hyperloglog & sketch) :
        sketch{std::addressof(sketch)},
        type{data_type::sketch}
    {}

    insert_iterator & operator=(uint64_t const value) noexcept
    {
        switch (type)
        {
        case data_type::unordered_set:
            assert(set != nullptr);
            set->emplace(value);
            break;
        case data_type::vector:
            assert(vec != nullptr);
            vec->emplace_back(value);
            break;
        case data_type::sketch:
            assert(sketch != nullptr);
            sketch->add(value);
            break;
        default:
#ifndef NDEBUG
            assert(false);
#else
            __builtin_unreachable();
#endif
        }
        return *this;
    }

    [[nodiscard]] constexpr insert_iterator & operator*() noexcept
    {
        return *this;
    }

    constexpr insert_iterator & operator++() noexcept
    {
        return *this;
    }

    constexpr insert_iterator operator++(int) noexcept
    {
        return *this;
    }

private:
    robin_hood::unordered_flat_set<uint64_t> * set{nullptr};
    std::vector<uint64_t> * vec{nullptr};
    sketch::hyperloglog * sketch{nullptr};

    enum class data_type : uint8_t
    {
        unordered_set,
        vector,
        sketch
    };

    data_type type{};
};

} // namespace seqan::hibf
