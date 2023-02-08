// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <type_traits>

#include <hibf/platform.hpp>

namespace hibf
{

template <typename t>
constexpr bool add_enum_bitwise_operators = false;

template <typename t>
constexpr t operator&(t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) & static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator|(t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) | static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator^(t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) ^ static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator~(t lhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(~static_cast<std::underlying_type_t<t>>(lhs));
}

template <typename t>
constexpr t & operator&=(t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs & rhs;
    return lhs;
}

template <typename t>
constexpr t & operator|=(t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs | rhs;
    return lhs;
}

template <typename t>
constexpr t & operator^=(t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs ^ rhs;
    return lhs;
}

} // namespace hibf
