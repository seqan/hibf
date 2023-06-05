// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/build/build_data.hpp>
#include <hibf/detail/layout/layout.hpp>
#include <hibf/detail/timer.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

namespace hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     hibf::interleaved_bloom_filter & ibf,
                     timer<concurrent::yes> & fill_ibf_timer);

void insert_into_ibf(build_data const & data,
                     layout::layout::user_bin const & record,
                     hibf::interleaved_bloom_filter & ibf);

} // namespace hibf
