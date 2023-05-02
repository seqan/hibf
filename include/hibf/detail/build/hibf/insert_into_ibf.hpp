// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/build/hibf/chopper_pack_record.hpp>
#include <hibf/interleaved_bloom_filter.hpp>

namespace hibf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     hibf::interleaved_bloom_filter<> & ibf)
{
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};
    size_t counter{};

    for (auto kmer : kmers)
    {
        if (counter == chunk_size)
        {
            counter = 0;
            ++chunk_number;
        }

        ibf.emplace(kmer, hibf::bin_index{bin_index + chunk_number});

        ++counter;
    }
}

template <typename config_type>
void insert_into_ibf(build_data<config_type> & data,
                     chopper_pack_record const & record,
                     hibf::interleaved_bloom_filter<> & ibf)
{
    auto const bin_index = hibf::bin_index{static_cast<size_t>(record.bin_indices.back())};

    for (auto && hash_sequence : data.hibf_config.input[record.user_bin_index])
        for (auto hash : hash_sequence)
            ibf.emplace(hash, bin_index);
}

} // namespace hibf
