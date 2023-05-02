// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <lemon/list_graph.h> // Must be first include.

#include <hibf/detail/build/hibf/build_data.hpp>
#include <hibf/detail/build/hibf/chopper_pack_record.hpp>

namespace hibf
{

template <typename config_type>
void update_user_bins(build_data<config_type> & data,
                      std::vector<int64_t> & filename_indices,
                      chopper_pack_record const & record)
{
    size_t const idx = record.user_bin_index;

    std::string & user_bin_filenames = data.hibf->user_bins.filename_of_user_bin(idx);
    for (auto const & filename : record.filenames)
    {
        user_bin_filenames += filename;
        user_bin_filenames += ';';
    }
    assert(!user_bin_filenames.empty());
    user_bin_filenames.pop_back();

    std::fill_n(filename_indices.begin() + record.bin_indices.back(), record.number_of_bins.back(), idx);
}

} // namespace hibf
