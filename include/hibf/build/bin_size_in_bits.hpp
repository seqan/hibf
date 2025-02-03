// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::bin_size_in_bits.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <cmath>   // for log, ceil, exp
#include <cstddef> // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf::build
{

/*!\brief Contains parameters for bin_size_in_bits.
 * \ingroup hibf_build
 * \qualifier strong
 */
struct bin_size_parameters
{
    double fpr{};
    size_t hash_count{};
    size_t elements{};
};

/*!\brief Computes the bin size.
 * \ingroup hibf_build
 */
inline constexpr size_t bin_size_in_bits(bin_size_parameters const & params)
{
    assert(params.hash_count > 0);
    assert(params.fpr > 0.0);
    assert(params.fpr < 1.0);

    double const numerator{-static_cast<double>(params.elements * params.hash_count)};
    double const denominator{std::log(1 - std::exp(std::log(params.fpr) / params.hash_count))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

} // namespace seqan::hibf::build
