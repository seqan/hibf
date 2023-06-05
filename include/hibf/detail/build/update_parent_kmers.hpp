// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements hibf::update_parent_kmers.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */
#pragma once

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/timer.hpp>
#include <hibf/platform.hpp>

namespace hibf
{

inline void update_parent_kmers(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                robin_hood::unordered_flat_set<uint64_t> const & kmers,
                                timer<concurrent::yes> & merge_kmers_timer)
{
    timer<concurrent::no> local_merge_kmers_timer{};
    local_merge_kmers_timer.start();
    parent_kmers.insert(kmers.begin(), kmers.end());
    local_merge_kmers_timer.stop();
    merge_kmers_timer += local_merge_kmers_timer;
}

} // namespace hibf
