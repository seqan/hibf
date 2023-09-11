// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides literals for unit conversions.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <hibf/platform.hpp>

namespace seqan::hibf::test
{

inline namespace literals
{

//!\brief A string literal that interprets the given number as MiB and converts to bytes.
static constexpr size_t operator""_MiB(unsigned long long int const number) noexcept
{
    return number << 20;
}

//!\brief A string literal that interprets the given number as GiB and converts to bytes.
static constexpr size_t operator""_GiB(unsigned long long int const number) noexcept
{
    return number << 30;
}

} // namespace literals

} // namespace seqan::hibf::test
