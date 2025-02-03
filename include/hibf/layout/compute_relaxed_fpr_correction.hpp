// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef> // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf::layout
{

/*!\brief Contains parameters for compute_relaxed_fpr_correction.
 * \ingroup hibf_layout
 * \qualifier strong
 */
struct relaxed_fpr_correction_parameters
{
    double fpr{};
    double relaxed_fpr{};
    size_t hash_count{};
};

/*!\brief Precompute size correction factor for merged bins which are allowed to have a relaxed FPR.
 * \ingroup hibf_layout
 */
[[nodiscard]] double compute_relaxed_fpr_correction(relaxed_fpr_correction_parameters const & params);

} // namespace seqan::hibf::layout
