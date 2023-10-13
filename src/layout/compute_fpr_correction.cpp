// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert> // for assert
#include <cmath>   // for log1p, exp, log
#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/layout/compute_fpr_correction.hpp> // for compute_fpr_correction
#include <hibf/misc/next_multiple_of_64.hpp>      // for next_multiple_of_64

namespace seqan::hibf::layout
{

std::vector<double> compute_fpr_correction(fpr_correction_parameters const & params)
{
    assert(params.fpr > 0.0 && params.fpr <= 1.0);
    assert(params.hash_count > 0u);

    size_t const max_tb = next_multiple_of_64(params.t_max);

    std::vector<double> fpr_correction(max_tb + 1u, 0.0);
    fpr_correction[1] = 1.0;

    // std::log1p(arg) = std::log(1 + arg). More precise than std::log(1 + arg) if arg is close to zero.
    double const numerator = std::log1p(-std::exp(std::log(params.fpr) / params.hash_count));

    for (size_t split = 2u; split <= max_tb; ++split)
    {
        double const log_target_fpr = std::log1p(-std::exp(std::log1p(-params.fpr) / split));
        fpr_correction[split] = numerator / std::log1p(-std::exp(log_target_fpr / params.hash_count));
        assert(fpr_correction[split] >= 1.0);
    }

    return fpr_correction;
}

} // namespace seqan::hibf::layout
