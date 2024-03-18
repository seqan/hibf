// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>  // for assert
#include <concepts> // for unsigned_integral
#include <cstddef>  // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns, for unsigned integral operands, `dividend / divisor` ceiled to the next integer value.
 * \ingroup hibf
 */
template <std::unsigned_integral t1, std::unsigned_integral t2>
[[nodiscard]] inline constexpr size_t divide_and_ceil(t1 const dividend, t2 const divisor) noexcept
{
    assert(divisor > 0u);
    assert(std::numeric_limits<size_t>::max() - divisor + 1u >= dividend); // Overflow detection
    return (static_cast<size_t>(dividend) + (divisor - 1u)) / divisor;
}

} // namespace seqan::hibf
