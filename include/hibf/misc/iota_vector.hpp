// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::iota_vector.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cassert>  // for assert
#include <concepts> // for unsigned_integral
#include <cstddef>  // for size_t
#include <limits>   // for numeric_limits
#include <numeric>  // for iota
#include <vector>   // for vector

#include <hibf/platform.hpp> // for HIBF_CONSTEXPR_VECTOR

namespace seqan::hibf
{

/*!\brief Creates a vector of size `size` with values from 0 to `size - 1`.
 * \tparam value_t The value type of the vector. Defaults to `size_t`.
 * \param[in] size The size of the vector.
 * \returns A vector of size `size` with values from 0 to `size - 1`.
*/
template <std::unsigned_integral value_t = size_t>
HIBF_CONSTEXPR_VECTOR std::vector<value_t> iota_vector(size_t const size)
{
    assert(size <= std::numeric_limits<value_t>::max());
    std::vector<value_t> result(size);
    std::iota(result.begin(), result.end(), value_t{});
    return result;
}

} // namespace seqan::hibf
