// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cstddef>

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns the number of empty bins that should be created by a given fraction of the total number of bins.
 * \param[in] tmax     The total number of bins.
 * \param[in] fraction The fraction of the total number of bins that should be empty.
 * \ingroup hibf
 * \sa https://godbolt.org/z/cMjbM39vj
 */
[[nodiscard]] constexpr size_t empty_bins_by_fraction(size_t const tmax, double const fraction) noexcept
{
    return std::clamp<size_t>(tmax * fraction, 1, tmax - 2) - (fraction == 0.0);
}

} // namespace seqan::hibf
