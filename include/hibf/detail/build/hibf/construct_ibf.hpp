// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/contrib/robin_hood.hpp>

#include <lemon/list_graph.h> /// Must be first include.

#include <seqan3/search/dream_index/detail/build/hibf/build_data.hpp>
#include <seqan3/search/dream_index/detail/build/hibf/insert_into_ibf.hpp>

namespace seqan3::hibf
{

template <seqan3::data_layout data_layout_mode, typename config_type>
seqan3::interleaved_bloom_filter<> construct_ibf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                                                 robin_hood::unordered_flat_set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data<data_layout_mode, config_type> & data,
                                                 bool is_root)
{
    auto bin_size_in_bits = [&] (size_t const number_of_kmers_to_be_stored)
    {
        double const numerator{-static_cast<double>(number_of_kmers_to_be_stored * data.hibf_config.number_of_hash_functions)};
        double const denominator{std::log(1 - std::exp(std::log(data.hibf_config.maximum_false_positive_rate) / data.hibf_config.number_of_hash_functions))};
        double const result{std::ceil(numerator / denominator)};
        return result;
    };

    auto & node_data = data.node_map[node];

    double const total_number_of_kmers{static_cast<double>(kmers.size())};
    size_t const kmers_per_bin{static_cast<size_t>(std::ceil(total_number_of_kmers / number_of_bins))};
    double const bin_bits{static_cast<double>(bin_size_in_bits(kmers_per_bin))};
    seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_bloom_filter<> ibf{bin_count, bin_size, seqan3::hash_function_count{data.hibf_config.number_of_hash_functions}};

    insert_into_ibf(parent_kmers, kmers, number_of_bins, node_data.max_bin_index, ibf, is_root);

    return ibf;
}

} // namespace seqan3::hibf
