// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements hibf::bin_size_in_bits.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cmath>   // for log, ceil, exp
#include <cstddef> // for size_t

#include <hibf/platform.hpp>

namespace hibf
{

inline size_t bin_size_in_bits(size_t const number_of_kmers_to_be_stored,
                               size_t const number_of_hash_functions,
                               double const maximum_false_positive_rate)
{
    double const numerator{-static_cast<double>(number_of_kmers_to_be_stored * number_of_hash_functions)};
    double const denominator{std::log(1 - std::exp(std::log(maximum_false_positive_rate) / number_of_hash_functions))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

} // namespace hibf
