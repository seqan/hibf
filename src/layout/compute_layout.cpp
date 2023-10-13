// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __sort_fn, sort
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for identity, function
#include <iterator>   // for inserter
#include <sstream>    // for basic_stringstream, stringstream
#include <utility>    // for addressof
#include <vector>     // for vector

#include <hibf/config.hpp>                        // for config
#include <hibf/contrib/robin_hood.hpp>            // for unordered_flat_set
#include <hibf/layout/compute_fpr_correction.hpp> // for compute_fpr_correction
#include <hibf/layout/compute_layout.hpp>         // for compute_layout
#include <hibf/layout/data_store.hpp>             // for data_store
#include <hibf/layout/hierarchical_binning.hpp>   // for hierarchical_binning
#include <hibf/layout/layout.hpp>                 // for layout
#include <hibf/sketch/estimate_kmer_counts.hpp>   // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>            // for hyperloglog

namespace seqan::hibf::layout
{

layout compute_layout(config const & config,
                      std::vector<size_t> const & kmer_counts,
                      std::vector<sketch::hyperloglog> const & sketches)
{
    layout resulting_layout{};

    // The output streams facilitate writing the layout file in hierarchical structure.
    // seqan::hibf::execute currently writes the filled buffers to the output file.
    std::stringstream output_buffer;
    std::stringstream header_buffer;

    data_store store{.false_positive_rate = config.maximum_false_positive_rate,
                     .hibf_layout = &resulting_layout,
                     .kmer_counts = std::addressof(kmer_counts),
                     .sketches = std::addressof(sketches)};
    store.fpr_correction = compute_fpr_correction({.fpr = config.maximum_false_positive_rate,
                                                   .hash_count = config.number_of_hash_functions,
                                                   .t_max = config.tmax});

    store.hibf_layout->top_level_max_bin_id = seqan::hibf::layout::hierarchical_binning{store, config}.execute();

    // sort records ascending by the number of bin indices (corresponds to the IBF levels)
    // GCOVR_EXCL_START
    std::ranges::sort(store.hibf_layout->max_bins,
                      [](auto const & r, auto const & l)
                      {
                          return r.previous_TB_indices.size() < l.previous_TB_indices.size();
                      });
    // GCOVR_EXCL_STOP

    return *store.hibf_layout;
}

} // namespace seqan::hibf::layout
