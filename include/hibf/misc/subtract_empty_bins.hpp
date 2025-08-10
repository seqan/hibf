// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm> // for clamp
#include <cmath>     // for ceil
#include <cstddef>   // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns the number of technical bins available for use.
 * \param[in] tmax     The total number of bins.
 * \param[in] fraction The fraction of the total number of bins that should be empty.
 * \ingroup hibf
 * \sa https://godbolt.org/z/cMjbM39vj
 */
[[nodiscard]] constexpr size_t subtract_empty_bins(size_t const tmax, double const fraction) noexcept
{
    // There must be at least 2 technical bins available without empty bins.
    // Otherwise, there would only ever be one technical bin available.
    if (fraction == 0.0 || tmax <= 2u)
        return tmax;

    size_t const number_of_empty_bins = std::clamp<size_t>(std::ceil(tmax * fraction), 1, tmax - 2);
    return tmax - number_of_empty_bins;
}

} // namespace seqan::hibf
