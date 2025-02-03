// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::update_parent_kmers.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */
#pragma once

#include <cstdint>    // for uint64_t
#include <functional> // for equal_to

#include <hibf/contrib/robin_hood.hpp> // for unordered_flat_set, hash
#include <hibf/misc/timer.hpp>         // for serial_timer, concurrent_timer
#include <hibf/platform.hpp>

namespace seqan::hibf::build
{

/*!\brief Updates stored values of the parent IBF.
 * \ingroup hibf_build
 */
inline void update_parent_kmers(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                robin_hood::unordered_flat_set<uint64_t> const & kmers,
                                concurrent_timer & merge_kmers_timer)
{
    serial_timer local_merge_kmers_timer{};
    local_merge_kmers_timer.start();
    parent_kmers.insert(kmers.begin(), kmers.end());
    local_merge_kmers_timer.stop();
    merge_kmers_timer += local_merge_kmers_timer;
}

} // namespace seqan::hibf::build
