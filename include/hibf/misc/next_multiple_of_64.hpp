// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef> // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns the smallest integer that is greater or equal to `value` and a multiple of 64. Returns 0 for value 0.
 * \param[in] value The Input value that is smaller or equal to the return value.
 * \ingroup hibf
 */
[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value) noexcept
{
    return ((value + 63) >> 6) << 6;
}

} // namespace seqan::hibf
