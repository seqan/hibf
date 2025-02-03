// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for shuffle, __fn, none_of
#include <cassert>    // for assert
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t
#include <functional> // for equal_to
#include <mutex>      // for mutex, lock_guard
#include <numeric>    // for iota
#include <optional>   // for optional
#include <random>     // for random_device, mt19937_64
#include <utility>    // for move
#include <vector>     // for vector, erase

#include <hibf/build/build_data.hpp>                      // for build_data
#include <hibf/build/compute_kmers.hpp>                   // for compute_kmers
#include <hibf/build/construct_ibf.hpp>                   // for construct_ibf
#include <hibf/build/insert_into_ibf.hpp>                 // for insert_into_ibf
#include <hibf/build/update_parent_kmers.hpp>             // for update_parent_kmers
#include <hibf/build/update_user_bins.hpp>                // for update_user_bins
#include <hibf/config.hpp>                                // for config
#include <hibf/contrib/robin_hood.hpp>                    // for unordered_flat_set, hash
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter, merged
#include <hibf/interleaved_bloom_filter.hpp>              // for interleaved_bloom_filter
#include <hibf/layout/compute_fpr_correction.hpp>         // for compute_fpr_correction
#include <hibf/layout/compute_layout.hpp>                 // for compute_layout
#include <hibf/layout/graph.hpp>                          // for graph
#include <hibf/layout/layout.hpp>                         // for layout
#include <hibf/misc/divide_and_ceil.hpp>                  // for divide_and_ceil
#include <hibf/misc/iota_vector.hpp>                      // for iota_vector
#include <hibf/misc/timer.hpp>                            // for concurrent_timer
#include <hibf/sketch/compute_sketches.hpp>               // for compute_sketches
#include <hibf/sketch/estimate_kmer_counts.hpp>           // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>                    // for hyperloglog

namespace seqan::hibf
{

size_t hierarchical_build(hierarchical_interleaved_bloom_filter & hibf,
                          robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                          layout::graph::node const & current_node,
                          build::build_data & data,
                          bool is_root,
                          size_t const parent_ibf_idx = 0u)
{
    size_t const ibf_pos{data.request_ibf_idx()};

    auto & technical_bin_to_ibf_id = hibf.next_ibf_id[ibf_pos];
    assert(technical_bin_to_ibf_id.empty());
    technical_bin_to_ibf_id.resize(current_node.number_of_technical_bins, ibf_pos);

    auto & technical_bin_to_user_bin_id = hibf.ibf_bin_to_user_bin_id[ibf_pos];
    assert(technical_bin_to_user_bin_id.empty());
    technical_bin_to_user_bin_id.resize(current_node.number_of_technical_bins, bin_kind::merged);

    auto & ibf = hibf.ibf_vector[ibf_pos];

    robin_hood::unordered_flat_set<uint64_t> kmers{};

    auto initialise_max_bin_kmers = [&]() -> size_t
    {
        if (current_node.max_bin_is_merged())
        {
            // recursively initialize favourite child first
            technical_bin_to_ibf_id[current_node.max_bin_index] =
                hierarchical_build(hibf,
                                   kmers,
                                   current_node.children[current_node.favourite_child_idx.value()],
                                   data,
                                   false,
                                   ibf_pos);
            return 1;
        }
        else // max bin is not a merged bin
        {
            // we assume that the max record is at the beginning of the list of remaining records.
            auto const & record = current_node.remaining_records[0];
            build::compute_kmers(kmers, data, record);
            build::update_user_bins(technical_bin_to_user_bin_id, record);

            return record.number_of_technical_bins;
        }
    };

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers();
    ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    auto loop_over_children = [&]()
    {
        if (current_node.children.empty())
            return;

        std::vector<layout::graph::node> children = current_node.children; // copy for threads

        size_t const number_of_mutex = divide_and_ceil(current_node.number_of_technical_bins, 64u);
        std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

        size_t number_of_threads{};
        std::vector<size_t> indices(children.size());
        std::iota(indices.begin(), indices.end(), size_t{});

        // We do not want to process the favourite child. It has already been processed prior.
        // https://godbolt.org/z/6Yav7hrG1
        if (current_node.max_bin_is_merged())
            std::erase(indices, current_node.favourite_child_idx.value());

        if (is_root)
        {
            // Shuffle indices: More likely to not block each other. Optimal: Interleave
            std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
            number_of_threads = data.config.threads;
        }
        else
        {
            number_of_threads = 1u;
        }

#pragma omp parallel for schedule(dynamic, 1) num_threads(number_of_threads)
        for (auto && index : indices)
        {
            auto & child = children[index];

            robin_hood::unordered_flat_set<uint64_t> local_kmers{};
            size_t const local_ibf_pos = hierarchical_build(hibf, local_kmers, child, data, false, ibf_pos);
            auto parent_bin_index = child.parent_bin_index;
            hibf.prev_ibf_id[local_ibf_pos] = {.ibf_idx = parent_ibf_idx, .bin_idx = parent_bin_index};
            {
                size_t const mutex_id{parent_bin_index / 64};
                std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                technical_bin_to_ibf_id[parent_bin_index] = local_ibf_pos;
                build::insert_into_ibf(local_kmers, 1, parent_bin_index, ibf, data.fill_ibf_timer);
                if (!is_root)
                    build::update_parent_kmers(parent_kmers, local_kmers, data.merge_kmers_timer);
            }
        }
    };

    loop_over_children();

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node.max_bin_is_merged()) ? 0u : 1u};
    for (size_t i = start; i < current_node.remaining_records.size(); ++i)
    {
        auto const & record = current_node.remaining_records[i];

        if (is_root && record.number_of_technical_bins == 1) // no splitting needed
        {
            build::insert_into_ibf(data, record, ibf);
        }
        else
        {
            compute_kmers(kmers, data, record);
            build::insert_into_ibf(kmers,
                                   record.number_of_technical_bins,
                                   record.storage_TB_id,
                                   ibf,
                                   data.fill_ibf_timer);
            if (!is_root)
                build::update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
        }

        build::update_user_bins(technical_bin_to_user_bin_id, record);
        kmers.clear();
    }

    return ibf_pos;
}

size_t hierarchical_build(hierarchical_interleaved_bloom_filter & hibf,
                          layout::graph::node const & root_node,
                          build::build_data & data)
{
    robin_hood::unordered_flat_set<uint64_t> root_kmers{};
    return hierarchical_build(hibf, root_kmers, root_node, data, true);
}

void build_index(hierarchical_interleaved_bloom_filter & hibf,
                 config const & config,
                 seqan::hibf::layout::layout const & hibf_layout)
{
    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1;

    hibf.ibf_vector.resize(number_of_ibfs);
    hibf.ibf_bin_to_user_bin_id.resize(number_of_ibfs);
    hibf.prev_ibf_id.resize(number_of_ibfs);
    hibf.next_ibf_id.resize(number_of_ibfs);

    build::build_data data{.config = config, .ibf_graph = {hibf_layout}};

    layout::graph::node const & root_node = data.ibf_graph.root;

    size_t const t_max{root_node.number_of_technical_bins};
    data.fpr_correction = layout::compute_fpr_correction({.fpr = config.maximum_fpr, //
                                                          .hash_count = config.number_of_hash_functions,
                                                          .t_max = t_max});

    hierarchical_build(hibf, root_node, data);

    // NOLINTBEGIN(performance-move-const-arg)
    hibf.index_allocation_timer = std::move(data.index_allocation_timer);
    hibf.user_bin_io_timer = std::move(data.user_bin_io_timer);
    hibf.merge_kmers_timer = std::move(data.merge_kmers_timer);
    hibf.fill_ibf_timer = std::move(data.fill_ibf_timer);
    // NOLINTEND(performance-move-const-arg)
}

hierarchical_interleaved_bloom_filter::hierarchical_interleaved_bloom_filter(config & configuration)
{
    configuration.validate_and_set_defaults();

    std::vector<sketch::hyperloglog> sketches{};
    std::vector<size_t> kmer_counts{};

    layout_compute_sketches_timer.start();
    sketch::compute_sketches(configuration, sketches);
    hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);
    layout_compute_sketches_timer.stop();

    // If rearrangement is enabled, i.e. seqan::hibf::config::disable_rearrangement is false:
    // `min_id == none` in seqan::hibf::sketch::toolbox::cluster_bins -> std::out_of_range "key not found"
    // Otherwise:
    // seqan::hibf::interleaved_bloom_filter constructor -> std::logic_error "The size of a bin must be > 0."
    assert(std::ranges::none_of(kmer_counts,
                                [](size_t const count)
                                {
                                    return count == 0u;
                                }));

    layout_dp_algorithm_timer.start();
    auto layout = layout::compute_layout(configuration,
                                         kmer_counts,
                                         sketches,
                                         iota_vector(configuration.number_of_user_bins),
                                         layout_union_estimation_timer,
                                         layout_rearrangement_timer);
    layout_dp_algorithm_timer.stop();

    number_of_user_bins = configuration.number_of_user_bins;
    build_index(*this, configuration, layout);
}

// seqan::hibf::config, seqan::hibf::layout::layout
hierarchical_interleaved_bloom_filter::hierarchical_interleaved_bloom_filter(config & configuration,
                                                                             layout::layout const & layout)
{
    configuration.validate_and_set_defaults();
    number_of_user_bins = configuration.number_of_user_bins;
    build_index(*this, configuration, layout);
}

} // namespace seqan::hibf
