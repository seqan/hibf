// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements seqan::hibf::compute_kmers.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <cstdint>    // for uint64_t
#include <functional> // for function

#include <hibf/build/build_data.hpp>    // for build_data
#include <hibf/build/compute_kmers.hpp> // for compute_kmers
#include <hibf/config.hpp>              // for insert_iterator, config
#include <hibf/contrib/robin_hood.hpp>  // for unordered_flat_set
#include <hibf/layout/layout.hpp>       // for layout
#include <hibf/misc/timer.hpp>          // for serial_timer, concurrent_timer

namespace seqan::hibf::build
{

void compute_kmers(robin_hood::unordered_flat_set<uint64_t> & kmers,
                   build_data const & data,
                   layout::layout::user_bin const & record)
{
    serial_timer local_user_bin_io_timer{};
    local_user_bin_io_timer.start();
    data.config.input_fn(record.idx, insert_iterator{kmers});
    local_user_bin_io_timer.stop();
    data.user_bin_io_timer += local_user_bin_io_timer;
}

} // namespace seqan::hibf::build
