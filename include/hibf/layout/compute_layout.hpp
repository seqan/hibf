#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/config.hpp>             // for config
#include <hibf/layout/layout.hpp>      // for layout
#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

namespace seqan::hibf::layout
{

/*!\brief Computes the layout.
 * \ingroup hibf_layout
 * \param[in] config The configuration to compute the layout with.
 * \param[in] kmer_counts The vector that will store the kmer counts (estimations).
 * \param[in] sketches The vector that will store the sketches.
 * \returns layout
 */
layout
compute_layout(config const & config, std::vector<size_t> const & kmer_counts, std::vector<sketch::hyperloglog> const & sketches);

} // namespace seqan::hibf::layout
