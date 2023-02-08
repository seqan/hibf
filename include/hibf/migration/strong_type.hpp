// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <concepts>
#include <type_traits>

#include <hibf/migration/add_enum_bitwise_operators.hpp>
// #include <seqan3/utility/type_traits/basic.hpp>

namespace hibf::detail
{

enum struct strong_type_skill
{
    none = 0,
    add = 1 << 0,
    subtract = 1 << 1,
    multiply = 1 << 2,
    divide = 1 << 3,
    modulo = 1 << 4,
    bitwise_and = 1 << 5,
    bitwise_or = 1 << 6,
    bitwise_xor = 1 << 7,
    bitwise_not = 1 << 8,
    bitwise_lshift = 1 << 9,
    bitwise_rshift = 1 << 10,
    logical_and = 1 << 11,
    logical_or = 1 << 12,
    logical_not = 1 << 13,
    increment = 1 << 14,
    decrement = 1 << 15,
    convert = 1 << 16,
    comparable = 1 << 17,
    additive = add | subtract,
    multiplicative = multiply | divide | modulo,
    bitwise_logic = bitwise_and | bitwise_or | bitwise_xor | bitwise_not,
    bitwise_shift = bitwise_lshift | bitwise_rshift,
    logic = logical_and | logical_or | logical_not
};

} //namespace hibf::detail

namespace hibf
{

template <>
constexpr bool add_enum_bitwise_operators<hibf::detail::strong_type_skill> = true;

} // namespace hibf

namespace hibf::detail
{

template <typename, typename, strong_type_skill>
class strong_type;

template <typename strong_type_t>
concept derived_from_strong_type =
    requires (strong_type_t && obj) {
        typename std::remove_reference_t<strong_type_t>::value_type;

        {
            std::remove_reference_t<strong_type_t>::skills
        };

        requires std::same_as<decltype(std::remove_reference_t<strong_type_t>::skills), strong_type_skill const>;

        requires std::derived_from<std::remove_cvref_t<strong_type_t>,
                                   strong_type<typename std::remove_reference_t<strong_type_t>::value_type,
                                               std::remove_cvref_t<strong_type_t>,
                                               std::remove_reference_t<strong_type_t>::skills>>;
    };

template <typename value_t, typename derived_t, strong_type_skill skills_ = strong_type_skill::none>
class strong_type
{
public:
    static constexpr strong_type_skill skills = skills_;
    using value_type = value_t;

    constexpr strong_type() noexcept = default;
    constexpr strong_type(strong_type const &) noexcept = default;
    constexpr strong_type(strong_type &&) noexcept = default;
    constexpr strong_type & operator=(strong_type const &) noexcept = default;
    constexpr strong_type & operator=(strong_type &&) noexcept = default;
    ~strong_type() noexcept = default;

    constexpr explicit strong_type(value_t _value) : value(std::move(_value))
    {}

    constexpr value_t & get() & noexcept
    {
        return value;
    }

    constexpr value_t const & get() const & noexcept
    {
        return value;
    }

    constexpr value_t && get() && noexcept
    {
        return std::move(value);
    }

    constexpr value_t const && get() const && noexcept
    {
        return std::move(value);
    }

    constexpr derived_t operator+(strong_type const & other)
        requires ((skills & strong_type_skill::add) != strong_type_skill::none)
    {
        return derived_t{get() + other.get()};
    }

    constexpr derived_t operator-(strong_type const & other)
        requires ((skills & strong_type_skill::subtract) != strong_type_skill::none)
    {
        return derived_t{get() - other.get()};
    }

    constexpr derived_t operator*(strong_type const & other)
        requires ((skills & strong_type_skill::multiply) != strong_type_skill::none)
    {
        return derived_t{get() * other.get()};
    }

    constexpr derived_t operator/(strong_type const & other)
        requires ((skills & strong_type_skill::divide) != strong_type_skill::none)
    {
        return derived_t{get() / other.get()};
    }

    constexpr derived_t operator%(strong_type const & other)
        requires ((skills & strong_type_skill::modulo) != strong_type_skill::none)
    {
        return derived_t{get() % other.get()};
    }

    constexpr derived_t operator&(strong_type const & other)
        requires ((skills & strong_type_skill::bitwise_and) != strong_type_skill::none)
    {
        return derived_t{get() & other.get()};
    }

    constexpr derived_t operator|(strong_type const & other)
        requires ((skills & strong_type_skill::bitwise_or) != strong_type_skill::none)
    {
        return derived_t{get() | other.get()};
    }

    constexpr derived_t operator^(strong_type const & other)
        requires ((skills & strong_type_skill::bitwise_xor) != strong_type_skill::none)
    {
        return derived_t{get() ^ other.get()};
    }

    constexpr derived_t operator~()
        requires ((skills & strong_type_skill::bitwise_not) != strong_type_skill::none)
    {
        return derived_t{~get()};
    }

    constexpr derived_t operator<<(strong_type const & other)
        requires ((skills & strong_type_skill::bitwise_lshift) != strong_type_skill::none)
    {
        return derived_t{get() << other.get()};
    }

    template <std::integral integral_t>
    constexpr derived_t operator<<(integral_t const shift)
        requires ((skills & strong_type_skill::bitwise_lshift) != strong_type_skill::none)
    {
        return derived_t{get() << shift};
    }

    constexpr derived_t operator>>(strong_type const & other)
        requires ((skills & strong_type_skill::bitwise_rshift) != strong_type_skill::none)
    {
        return derived_t{get() >> other.get()};
    }

    template <std::integral integral_t>
    constexpr derived_t operator>>(integral_t const shift)
        requires ((skills & strong_type_skill::bitwise_rshift) != strong_type_skill::none)
    {
        return derived_t{get() >> shift};
    }

    constexpr bool operator&&(strong_type const & other)
        requires ((skills & strong_type_skill::logical_and) != strong_type_skill::none)
    {
        return get() && other.get();
    }

    constexpr bool operator||(strong_type const & other)
        requires ((skills & strong_type_skill::logical_or) != strong_type_skill::none)
    {
        return get() || other.get();
    }

    constexpr bool operator!()
        requires ((skills & strong_type_skill::logical_not) != strong_type_skill::none)
    {
        return !get();
    }

    constexpr derived_t & operator++()
        requires ((skills & strong_type_skill::increment) != strong_type_skill::none)
    {
        ++get();
        return static_cast<derived_t &>(*this);
    }

    constexpr derived_t operator++(int)
        requires ((skills & strong_type_skill::increment) != strong_type_skill::none)
    {
        derived_t tmp{get()};
        ++get();
        return tmp;
    }

    constexpr derived_t & operator--()
        requires ((skills & strong_type_skill::decrement) != strong_type_skill::none)
    {
        --get();
        return static_cast<derived_t &>(*this);
    }

    constexpr derived_t operator--(int)
        requires ((skills & strong_type_skill::decrement) != strong_type_skill::none)
    {
        derived_t tmp{get()};
        --get();
        return tmp;
    }

    constexpr bool operator==(strong_type const & rhs) const
        requires ((skills & strong_type_skill::comparable) != strong_type_skill::none)
    {
        return get() == rhs.get();
    }

    constexpr bool operator!=(strong_type const & rhs) const
        requires ((skills & strong_type_skill::comparable) != strong_type_skill::none)
    {
        return !(*this == rhs);
    }

    explicit constexpr operator value_t() const
        requires ((skills & strong_type_skill::convert) != strong_type_skill::none)
    {
        return get();
    }

private:
    value_t value;
};

} // namespace hibf::detail
