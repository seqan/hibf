// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements hibf::compute_kmers.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/detail/build/compute_kmers.hpp>

namespace hibf
{

void compute_kmers(robin_hood::unordered_flat_set<uint64_t> & kmers,
                   build_data const & data,
                   layout::layout::user_bin const & record)
{
    timer<concurrent::no> local_user_bin_io_timer{};
    local_user_bin_io_timer.start();
    data.hibf_config.input_fn(record.idx, std::inserter(kmers, kmers.begin()));
    local_user_bin_io_timer.stop();
    data.user_bin_io_timer += local_user_bin_io_timer;
}

} // namespace hibf