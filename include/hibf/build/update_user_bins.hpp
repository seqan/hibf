// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm> // for fill_n
#include <cstdint>   // for uint64_t
#include <vector>    // for vector

#include <hibf/layout/layout.hpp> // for layout

namespace seqan::hibf::build
{

/*!\brief Updates user bins stored in HIBF.
 * \ingroup hibf_build
 */
inline void update_user_bins(std::vector<uint64_t> & technical_bin_to_user_bin_id,
                             layout::layout::user_bin const & record)
{
    std::fill_n(technical_bin_to_user_bin_id.begin() + record.storage_TB_id,
                record.number_of_technical_bins,
                record.idx);
}

} // namespace seqan::hibf::build
