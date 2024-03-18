// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/platform.hpp>

namespace seqan::hibf::layout
{

/*!\brief Contains parameters for compute_fpr_correction.
 * \ingroup hibf_layout
 * \qualifier strong
 */
struct fpr_correction_parameters
{
    double fpr{};
    size_t hash_count{};
    size_t t_max{};
};

/*!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
 * \ingroup hibf_layout
 * \sa https://godbolt.org/z/zTj1v9W94
 */
std::vector<double> compute_fpr_correction(fpr_correction_parameters const & params);

} // namespace seqan::hibf::layout
