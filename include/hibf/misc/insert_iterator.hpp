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
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/platform.hpp>
#include <hibf/sketch/hyperloglog.hpp>

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

    constexpr insert_iterator() = default;
    constexpr insert_iterator(insert_iterator const &) = default;
    constexpr insert_iterator(insert_iterator &&) = default;
    constexpr insert_iterator & operator=(insert_iterator const &) = default;
    constexpr insert_iterator & operator=(insert_iterator &&) = default;
    constexpr ~insert_iterator() = default;

    using set_t = robin_hood::unordered_flat_set<uint64_t>;
    using sketch_t = sketch::hyperloglog;
    using ibf_t = interleaved_bloom_filter;
    using function_t = std::function<void(uint64_t const)>;

    explicit constexpr insert_iterator(set_t & set) : ptr{std::addressof(set)}, type{data_type::unordered_set}
    {}

    explicit constexpr insert_iterator(sketch_t & sketch) : ptr{std::addressof(sketch)}, type{data_type::sketch}
    {}

    explicit constexpr insert_iterator(ibf_t & ibf, size_t ibf_bin_index) :
        ptr{std::addressof(ibf)},
        ibf_bin_index{ibf_bin_index},
        type{data_type::ibf}
    {}

    explicit constexpr insert_iterator(function_t & fun) : ptr{std::addressof(fun)}, type{data_type::function}
    {}

    [[gnu::always_inline, gnu::flatten]] inline insert_iterator & operator=(uint64_t const value) noexcept
    {
        assert(ptr != nullptr);

        switch (type)
        {
        case data_type::unordered_set:
            static_cast<set_t *>(ptr)->emplace(value);
            break;
        case data_type::sketch:
            static_cast<sketch_t *>(ptr)->add(value);
            break;
        case data_type::ibf:
            static_cast<ibf_t *>(ptr)->emplace(value, static_cast<bin_index>(ibf_bin_index));
            break;
        default:
            assert(type == data_type::function);
            static_cast<function_t *>(ptr)->operator()(value);
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
    void * ptr{nullptr};

    enum class data_type : uint8_t
    {
        unordered_set,
        sketch,
        ibf,
        function
    };

    size_t ibf_bin_index{};
    data_type type{};
};

} // namespace seqan::hibf
