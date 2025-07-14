// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>   // for floor
#include <cstddef> // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns the total number of bins such that a `fraction` of bins is empty, and `tmax` many bins are non-empty.
 * \param[in] tmax     The number of of non-empty bins.
 * \param[in] fraction The fraction of the total number of bins that should be empty.
 * \ingroup hibf
 */
[[nodiscard]] constexpr size_t add_empty_bins(size_t const tmax, double const fraction) noexcept
{
    return std::floor(tmax / (1.0 - fraction));
}

} // namespace seqan::hibf
