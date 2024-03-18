// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>   // for assert
#include <cinttypes> // for uint64_t
#include <cstddef>   // for ptrdiff_t
#include <iterator>  // for output_iterator_tag
#include <memory>    // for addressof
#include <vector>    // for vector

#include <hibf/contrib/robin_hood.hpp> // for unordered_flat_set
#include <hibf/platform.hpp>

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
        is_set{true}
    {}

    explicit constexpr insert_iterator(std::vector<uint64_t> & vec) : vec{std::addressof(vec)}, is_set{false}
    {}

    insert_iterator & operator=(uint64_t const value) noexcept
    {
        if (is_set)
        {
            assert(set != nullptr);
            set->emplace(value);
        }
        else
        {
            assert(vec != nullptr);
            vec->emplace_back(value);
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
    bool is_set{false};
};

} // namespace seqan::hibf
