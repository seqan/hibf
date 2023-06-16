#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/config.hpp>                    // for config
#include <hibf/detail/layout/layout.hpp>      // for layout
#include <hibf/detail/sketch/hyperloglog.hpp> // for hyperloglog

namespace hibf::layout
{

/*!\brief Computes the layout and stores the kmer_counts and sketches in the respective vectors for further use.
 * \param hibf_config The configuration to compute the layout with.
 * \param[in,out] kmer_counts The vector that will store the kmer counts (estimations).
 * \param[in,out] sketches The vector that will store the sketches.
 * \returns layout
 */
layout compute_layout(config const & hibf_config,
                      std::vector<size_t> & kmer_counts,
                      std::vector<sketch::hyperloglog> & sketches);

layout compute_layout(config const & hibf_config);

} // namespace hibf::layout