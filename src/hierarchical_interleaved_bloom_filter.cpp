// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>  // for fill_n, max, shuffle
#include <cinttypes>  // for uint64_t, int64_t
#include <cmath>      // for ceil, sqrt
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iostream>   // for char_traits, operator<<, basic_ostream, cerr
#include <mutex>      // for mutex, lock_guard
#include <numeric>    // for iota
#include <optional>   // for optional
#include <random>     // for random_device, mt19937_64
#include <stdexcept>  // for invalid_argument
#include <utility>    // for move
#include <vector>     // for vector, erase

#include <hibf/build/build_data.hpp>                      // for build_data
#include <hibf/build/compute_kmers.hpp>                   // for compute_kmers
#include <hibf/build/construct_ibf.hpp>                   // for construct_ibf
#include <hibf/build/insert_into_ibf.hpp>                 // for insert_into_ibf
#include <hibf/build/update_parent_kmers.hpp>             // for update_parent_kmers
#include <hibf/build/update_user_bins.hpp>                // for update_user_bins
#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/contrib/robin_hood.hpp>                    // for unordered_flat_set
#include <hibf/detail/timer.hpp>                          // for timer
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/interleaved_bloom_filter.hpp>              // for interleaved_bloom_filter
#include <hibf/layout/compute_fpr_correction.hpp>         // for compute_fpr_correction
#include <hibf/layout/compute_layout.hpp>                 // for compute_layout
#include <hibf/layout/graph.hpp>                          // for graph
#include <hibf/layout/layout.hpp>                         // for layout
#include <hibf/next_multiple_of_64.hpp>                   // for next_multiple_of_64
#include <hibf/user_bins_type.hpp>                        // for user_bins_type

namespace seqan::hibf
{

size_t hierarchical_build(hierarchical_interleaved_bloom_filter & hibf,
                          robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                          layout::graph::node const & current_node,
                          build::build_data & data,
                          bool is_root)
{
    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node.number_of_technical_bins, -1);
    robin_hood::unordered_flat_set<uint64_t> kmers{};

    auto initialise_max_bin_kmers = [&]() -> size_t
    {
        if (current_node.favourite_child_idx.has_value()) // max bin is a merged bin
        {
            // recursively initialize favourite child first
            ibf_positions[current_node.max_bin_index] =
                hierarchical_build(hibf,
                                   kmers,
                                   current_node.children[current_node.favourite_child_idx.value()],
                                   data,
                                   false);
            return 1;
        }
        else // max bin is not a merged bin
        {
            // we assume that the max record is at the beginning of the list of remaining records.
            auto const & record = current_node.remaining_records[0];
            build::compute_kmers(kmers, data, record);
            build::update_user_bins(filename_indices, record);

            return record.number_of_technical_bins;
        }
    };

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers();
    auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    auto loop_over_children = [&]()
    {
        if (current_node.children.empty())
            return;

        std::vector<layout::graph::node> children = current_node.children; // copy for threads

        size_t const number_of_mutex = (current_node.number_of_technical_bins + 63) / 64;
        std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

        size_t number_of_threads{};
        std::vector<size_t> indices(children.size());
        std::iota(indices.begin(), indices.end(), size_t{});

        // We do not want to process the favourite child. It has already been processed prior.
        // https://godbolt.org/z/6Yav7hrG1
        if (current_node.favourite_child_idx.has_value())
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

            robin_hood::unordered_flat_set<uint64_t> kmers{};
            size_t const ibf_pos = hierarchical_build(hibf, kmers, child, data, false);
            auto parent_bin_index = child.parent_bin_index;
            {
                size_t const mutex_id{parent_bin_index / 64};
                std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                ibf_positions[parent_bin_index] = ibf_pos;
                build::insert_into_ibf(kmers, 1, parent_bin_index, ibf, data.fill_ibf_timer);
                if (!is_root)
                    build::update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
            }
        }
    };

    loop_over_children();

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node.favourite_child_idx.has_value()) ? 0u : 1u};
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

        build::update_user_bins(filename_indices, record);
        kmers.clear();
    }

    hibf.ibf_vector[ibf_pos] = std::move(ibf);
    hibf.next_ibf_id[ibf_pos] = std::move(ibf_positions);
    hibf.user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

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
                 seqan::hibf::layout::layout & hibf_layout)
{
    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1;

    hibf.ibf_vector.resize(number_of_ibfs);
    hibf.user_bins.set_ibf_count(number_of_ibfs);
    hibf.user_bins.set_user_bin_count(hibf_layout.user_bins.size());
    hibf.next_ibf_id.resize(number_of_ibfs);

    build::build_data data{.config = config, .ibf_graph = {hibf_layout}};

    layout::graph::node const & root_node = data.ibf_graph.root;

    size_t const t_max{root_node.number_of_technical_bins};
    data.fpr_correction = layout::compute_fpr_correction(
        {.fpr = config.maximum_false_positive_rate, .hash_count = config.number_of_hash_functions, .t_max = t_max});

    hierarchical_build(hibf, root_node, data);

    hibf.index_allocation_timer = std::move(data.index_allocation_timer);
    hibf.user_bin_io_timer = std::move(data.user_bin_io_timer);
    hibf.merge_kmers_timer = std::move(data.merge_kmers_timer);
    hibf.fill_ibf_timer = std::move(data.fill_ibf_timer);
}

/*!\brief Checks several variables of seqan::hibf::config and sets default values if necessary.
 *
 * The following checks are performed and might throw an exception if the configuration doesn't pass validation:
 * * If seqan::hibf::config::number_of_user_bins is `0` it is considered `not set` and an exception will be thrown
 *   because this paprameter is required.
 * * If seqan::hibf::config::tmax is `0` and seqan::hibf::config::number_of_user_bins is `>= 1ULL << 32` an exception
 *   will be thrown because a default tmax cannot be computed.
 *
 * The configuration might be modified as follows before passed to the HIBF construction algorithm:
 * * If seqan::hibf::config::disable_estimate_union is set to true but seqan::hibf::config::disable_rearrangement
 *   is not, seqan::hibf::config::disable_rearrangement will be fordced to be true also. Without union estimations,
 *   no rearrangement can be done.
 * * If seqan::hibf::config::tmax is `0` it is considered `not set` and a default will be computed via
 *   `static_cast<uint16_t>(std::ceil(std::sqrt(cfg.number_of_user_bins)))`.
 * * If seqan::hibf::config::tmax is **not** `0` but also not a multiple of 64, it is increased to the next multiple of
 *   64 to avoid uneccessary space consumption (e.g. value `60` will be increased to `64` or `1000` increased to `1024`).
 */
void check_config_and_set_defaults(config & cfg)
{
    if (cfg.disable_estimate_union)
        cfg.disable_rearrangement = true;

    if (cfg.number_of_user_bins = 0)
        throw std::invalid_argument{
            "[HIBF CONFIG ERROR] You didn't set config::number_of_user_bins but it's required."};

    if (cfg.tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (cfg.number_of_user_bins >= 1ULL << 32) // sqrt is bigger than uint16_t
        {
            throw std::invalid_argument{
                "[HIBF CONFIG ERROR] Too many user-bins/samples to compute a default tmax. " // GCOVR_EXCL_LINE
                "Please set a tmax manually."};                                              // GCOVR_EXCL_LINE
        }
        else
        {
            cfg.tmax =
                seqan::hibf::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(cfg.number_of_user_bins))));
        }
    }
    else if (cfg.tmax % 64 != 0)
    {
        cfg.tmax = seqan::hibf::next_multiple_of_64(cfg.tmax);
        std::cerr << "[HIBF CONFIG WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << cfg.tmax << ".\n";
    }
}

hierarchical_interleaved_bloom_filter::hierarchical_interleaved_bloom_filter(config & configuration)
{
    check_config_and_set_defaults(configuration);
    auto layout = layout::compute_layout(configuration);
    build_index(*this, configuration, layout);
}

hierarchical_interleaved_bloom_filter::hierarchical_interleaved_bloom_filter(
    std::function<void(size_t const, insert_iterator &&)> input_fn,
    std::istream & layout_stream,
    size_t const threads)
{
    // read config and layout from file
    config configuration;
    layout::layout hibf_layout;
    configuration.read_from(layout_stream);
    hibf_layout.read_from(layout_stream);

    configuration.input_fn = input_fn; // set input as it cannot be serialized.
    configuration.threads = threads;   // set threads as it shouldn't use what was used to compute the layout.

    build_index(*this, configuration, hibf_layout);
}

} // namespace seqan::hibf
