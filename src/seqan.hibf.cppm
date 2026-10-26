// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides the `seqan.hibf` C++20 module.
 *
 * \details
 *
 * All headers are included in the *global module fragment*. Consequently, every entity of the library stays attached
 * to the global module, exactly as it is in a header build:
 *
 * * `libhibf.a` is bit-identical whether or not the module is built; the module adds one (almost empty) object file.
 * * The source files in `src` need no changes: they are ordinary translation units, and the symbols they define are
 *   the same symbols the module declares.
 * * A translation unit may `#include <hibf/...>` and `import seqan.hibf;` at the same time. Both refer to the same
 *   entities, so there is no ODR hazard and downstream projects can migrate file by file.
 *
 * The price is that the exported names have to be listed explicitly below. Names that are absent from this list are
 * still *reachable* (e.g. a member function of an exported class can use them), just not *visible* to importers.
 *
 * \attention GCC (as of 16.2) rejects a textual `#include` of a standard library header that appears *after* an
 *            `import seqan.hibf;` in the same translation unit. Put includes before imports. Clang has no such
 *            restriction.
 */

module;

// ============================================================================
//  Global module fragment
// ============================================================================

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/build/build_data.hpp>
#include <hibf/build/compute_kmers.hpp>
#include <hibf/build/construct_ibf.hpp>
#include <hibf/build/insert_into_ibf.hpp>
#include <hibf/build/update_parent_kmers.hpp>
#include <hibf/build/update_user_bins.hpp>
#include <hibf/cereal/concepts.hpp>
#include <hibf/cereal/path.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/compute_relaxed_fpr_correction.hpp>
#include <hibf/layout/data_store.hpp>
#include <hibf/layout/graph.hpp>
#include <hibf/layout/hierarchical_binning.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/layout/prefixes.hpp>
#include <hibf/layout/print_matrix.hpp>
#include <hibf/layout/simple_binning.hpp>
#include <hibf/misc/add_empty_bins.hpp>
#include <hibf/misc/bit_vector.hpp>
#include <hibf/misc/counting_vector.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/misc/insert_iterator.hpp>
#include <hibf/misc/iota_vector.hpp>
#include <hibf/misc/next_multiple_of_64.hpp>
#include <hibf/misc/print.hpp>
#include <hibf/misc/subtract_empty_bins.hpp>
#include <hibf/misc/timer.hpp>
#include <hibf/misc/unreachable.hpp>
#include <hibf/platform.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>
#include <hibf/sketch/toolbox.hpp>
#include <hibf/sketch/wyhash.hpp>
#include <hibf/version.hpp>

// ============================================================================
//  Module purview
// ============================================================================

export module seqan.hibf;

//!\brief The main namespace of the HIBF library.
export namespace seqan::hibf
{

// version.hpp
using ::seqan::hibf::hibf_version;
using ::seqan::hibf::hibf_version_cstring;
using ::seqan::hibf::hibf_version_major;
using ::seqan::hibf::hibf_version_minor;
using ::seqan::hibf::hibf_version_patch;

// cereal/concepts.hpp
using ::seqan::hibf::cereal_archive;
using ::seqan::hibf::cereal_input_archive;
using ::seqan::hibf::cereal_output_archive;
using ::seqan::hibf::cereal_text_archive;

// config.hpp
using ::seqan::hibf::config;

// interleaved_bloom_filter.hpp
using ::seqan::hibf::bin_count;
using ::seqan::hibf::bin_index;
using ::seqan::hibf::bin_size;
using ::seqan::hibf::empty_bin_fraction;
using ::seqan::hibf::hash_function_count;
using ::seqan::hibf::inspector;
using ::seqan::hibf::interleaved_bloom_filter;
using ::seqan::hibf::track_occupancy;

// hierarchical_interleaved_bloom_filter.hpp
using ::seqan::hibf::hierarchical_interleaved_bloom_filter;

// misc/
using ::seqan::hibf::add_empty_bins;
using ::seqan::hibf::bit_vector;
using ::seqan::hibf::concurrent_timer;
using ::seqan::hibf::counting_vector;
using ::seqan::hibf::divide_and_ceil;
using ::seqan::hibf::insert_iterator;
using ::seqan::hibf::iota_vector;
using ::seqan::hibf::next_multiple_of_64;
using ::seqan::hibf::print;
using ::seqan::hibf::print_t;
using ::seqan::hibf::serial_timer;
using ::seqan::hibf::subtract_empty_bins;
using ::seqan::hibf::unreachable;

} // namespace seqan::hibf

//!\brief Constants that indicate the kind of a technical bin.
export namespace seqan::hibf::bin_kind
{

using ::seqan::hibf::bin_kind::deleted;
using ::seqan::hibf::bin_kind::merged;

} // namespace seqan::hibf::bin_kind

//!\brief The layout algorithms.
export namespace seqan::hibf::layout
{

using ::seqan::hibf::layout::compute_fpr_correction;
using ::seqan::hibf::layout::compute_layout;
using ::seqan::hibf::layout::compute_relaxed_fpr_correction;
using ::seqan::hibf::layout::data_store;
using ::seqan::hibf::layout::fpr_correction_parameters;
using ::seqan::hibf::layout::graph;
using ::seqan::hibf::layout::hierarchical_binning;
using ::seqan::hibf::layout::layout;
using ::seqan::hibf::layout::print_matrix;
using ::seqan::hibf::layout::relaxed_fpr_correction_parameters;
using ::seqan::hibf::layout::simple_binning;

} // namespace seqan::hibf::layout

//!\brief Prefixes used when writing a layout file.
export namespace seqan::hibf::prefix
{

using ::seqan::hibf::prefix::layout_column_names;
using ::seqan::hibf::prefix::layout_first_header_line;
using ::seqan::hibf::prefix::layout_fullest_technical_bin_idx;
using ::seqan::hibf::prefix::layout_header;
using ::seqan::hibf::prefix::layout_lower_level;
using ::seqan::hibf::prefix::layout_top_level;
using ::seqan::hibf::prefix::meta_header;
using ::seqan::hibf::prefix::meta_hibf_config_end;
using ::seqan::hibf::prefix::meta_hibf_config_start;

} // namespace seqan::hibf::prefix

//!\brief Building an HIBF from a layout.
export namespace seqan::hibf::build
{

using ::seqan::hibf::build::bin_size_in_bits;
using ::seqan::hibf::build::bin_size_parameters;
using ::seqan::hibf::build::build_data;
using ::seqan::hibf::build::compute_kmers;
using ::seqan::hibf::build::construct_ibf;
using ::seqan::hibf::build::insert_into_ibf;
using ::seqan::hibf::build::update_parent_kmers;
using ::seqan::hibf::build::update_user_bins;

} // namespace seqan::hibf::build

//!\brief Sketches used to estimate cardinalities for the layout algorithm.
export namespace seqan::hibf::sketch
{

using ::seqan::hibf::sketch::compute_sketches;
using ::seqan::hibf::sketch::estimate_kmer_counts;
using ::seqan::hibf::sketch::hyperloglog;
using ::seqan::hibf::sketch::minhashes;
using ::seqan::hibf::sketch::wyhash;

} // namespace seqan::hibf::sketch

//!\brief Union estimation and similarity rearrangement.
export namespace seqan::hibf::sketch::toolbox
{

using ::seqan::hibf::sketch::toolbox::cluster_bins;
using ::seqan::hibf::sketch::toolbox::clustering_node;
using ::seqan::hibf::sketch::toolbox::distance_matrix;
using ::seqan::hibf::sketch::toolbox::entry;
// `estimate_interval` is behind `#if 0` in hibf/sketch/toolbox.hpp and therefore not exported.
using ::seqan::hibf::sketch::toolbox::neighbor;
using ::seqan::hibf::sketch::toolbox::precompute_initial_union_estimates;
using ::seqan::hibf::sketch::toolbox::precompute_union_estimates_for;
using ::seqan::hibf::sketch::toolbox::prio_queue;
using ::seqan::hibf::sketch::toolbox::prune;
using ::seqan::hibf::sketch::toolbox::random_shuffle;
using ::seqan::hibf::sketch::toolbox::rearrange_bins;
using ::seqan::hibf::sketch::toolbox::rotate;
using ::seqan::hibf::sketch::toolbox::sort_by_cardinalities;
using ::seqan::hibf::sketch::toolbox::trace;

} // namespace seqan::hibf::sketch::toolbox

/*!\brief Serialisation support for `std::filesystem::path`, provided by hibf/cereal/path.hpp.
 * \details
 * cereal finds `save`/`load` via ADL on `cereal::` types. The overloads must therefore be visible to importers, not
 * merely reachable. Do not drop this block: Clang happens to compile without it, GCC does not.
 */
export namespace cereal
{

using ::cereal::CEREAL_LOAD_FUNCTION_NAME;
using ::cereal::CEREAL_SAVE_FUNCTION_NAME;

} // namespace cereal
