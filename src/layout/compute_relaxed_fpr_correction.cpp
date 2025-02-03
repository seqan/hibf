// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert> // for assert
#include <cmath>   // for log1p, exp, log

#include <hibf/layout/compute_relaxed_fpr_correction.hpp> // for compute_fpr_correction

namespace seqan::hibf::layout
{

double compute_relaxed_fpr_correction(relaxed_fpr_correction_parameters const & params)
{
    assert(params.fpr > 0.0 && params.fpr <= 1.0);
    assert(params.relaxed_fpr > 0.0 && params.relaxed_fpr <= 1.0);
    assert(params.hash_count > 0u);
    assert(params.fpr <= params.relaxed_fpr);

    double const numerator = std::log1p(-std::exp(std::log(params.fpr) / params.hash_count));
    double const denominator = std::log1p(-std::exp(std::log(params.relaxed_fpr) / params.hash_count));
    double const relaxed_fpr_correction = numerator / denominator;

    assert(relaxed_fpr_correction > 0.0 && relaxed_fpr_correction <= 1.0);
    return relaxed_fpr_correction;
}

} // namespace seqan::hibf::layout
