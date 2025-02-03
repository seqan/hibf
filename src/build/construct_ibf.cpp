// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert>    // for assert
#include <cmath>      // for ceil
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t
#include <functional> // for equal_to
#include <vector>     // for vector

#include <hibf/build/bin_size_in_bits.hpp>    // for bin_size_in_bits
#include <hibf/build/build_data.hpp>          // for build_data
#include <hibf/build/construct_ibf.hpp>       // for construct_ibf
#include <hibf/build/insert_into_ibf.hpp>     // for insert_into_ibf
#include <hibf/build/update_parent_kmers.hpp> // for update_parent_kmers
#include <hibf/config.hpp>                    // for config
#include <hibf/contrib/robin_hood.hpp>        // for unordered_flat_set, hash
#include <hibf/interleaved_bloom_filter.hpp>  // for interleaved_bloom_filter, bin_count, bin_size, hash_function_c...
#include <hibf/layout/graph.hpp>              // for graph
#include <hibf/misc/divide_and_ceil.hpp>      // for divide_and_ceil
#include <hibf/misc/timer.hpp>                // for serial_timer, concurrent_timer

namespace seqan::hibf::build
{

seqan::hibf::interleaved_bloom_filter construct_ibf(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                                    robin_hood::unordered_flat_set<uint64_t> & kmers,
                                                    size_t const number_of_bins,
                                                    layout::graph::node const & ibf_node,
                                                    build_data & data,
                                                    bool is_root)
{
    bool const max_bin_is_merged = ibf_node.max_bin_is_merged();
    assert(!max_bin_is_merged || number_of_bins == 1u); // merged max bin implies (=>) number of bins == 1

    size_t const kmers_per_bin = divide_and_ceil(kmers.size(), number_of_bins);
    double const fpr = max_bin_is_merged ? data.config.relaxed_fpr : data.config.maximum_fpr;

    size_t const bin_bits{bin_size_in_bits({.fpr = fpr, //
                                            .hash_count = data.config.number_of_hash_functions,
                                            .elements = kmers_per_bin})};
    // data.fpr_correction[1] == 1.0, but we can avoid floating point operations with the ternary.
    // Check number_of_bins instead of max_bin_is_merged, because split bins can also occupy only one technical bin.
    seqan::hibf::bin_size const bin_size{
        number_of_bins == 1u ? bin_bits
                             : static_cast<size_t>(std::ceil(bin_bits * data.fpr_correction[number_of_bins]))};
    seqan::hibf::bin_count const bin_count{ibf_node.number_of_technical_bins};

    serial_timer local_index_allocation_timer{};
    local_index_allocation_timer.start();
    seqan::hibf::interleaved_bloom_filter ibf{bin_count,
                                              bin_size,
                                              seqan::hibf::hash_function_count{data.config.number_of_hash_functions},
                                              data.config.empty_bin_fraction > 0.0};

    local_index_allocation_timer.stop();
    data.index_allocation_timer += local_index_allocation_timer;

    insert_into_ibf(kmers, number_of_bins, ibf_node.max_bin_index, ibf, data.fill_ibf_timer);
    if (!is_root)
        update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);

    return ibf;
}

} // namespace seqan::hibf::build
