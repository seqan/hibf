// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::compute_kmers.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cinttypes> // for uint64_t

#include <hibf/build/build_data.hpp>   // for build_data
#include <hibf/contrib/robin_hood.hpp> // for unordered_flat_set
#include <hibf/layout/layout.hpp>      // for layout

namespace seqan::hibf::build
{

/*!\brief Computes kmers.
 * \ingroup hibf_build
 */
void compute_kmers(robin_hood::unordered_flat_set<uint64_t> & kmers,
                   build_data const & data,
                   layout::layout::user_bin const & record);

} // namespace seqan::hibf::build
