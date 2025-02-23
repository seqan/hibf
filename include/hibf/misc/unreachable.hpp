// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::unreachable.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

// IWYU pragma: begin_exports

#include <cassert> // for assert
#include <utility> // for std::unreachable

// IWYU pragma: end_exports

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\fn inline void seqan::hibf::unreachable()
 * \brief Marks unreachable code paths.
 * \details
 * In debug mode, it triggers an assertion failure.
 * In release mode, it calls `std::unreachable`.
 * ### Example
 * \include test/snippet/unreachable.cpp
 */
#ifndef NDEBUG
inline void unreachable() // GCOVR_EXCL_LINE
{
    assert(false); // GCOVR_EXCL_LINE
}
#else
using std::unreachable;
#endif

} // namespace seqan::hibf
