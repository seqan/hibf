#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/config.hpp>             // for config
#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

namespace seqan::hibf::sketch
{

/*!\brief Computes the kmer_counts and sketches and stores them in the respective vectors for further use.
 * \ingroup hibf_layout
 * \param[in] config The configuration to compute the layout with.
 * \param[in,out] kmer_counts The vector that will store the kmer counts (estimations).
 * \param[in,out] sketches The vector that will store the sketches.
 */
void compute_sketches(config const & config,
                      std::vector<size_t> & kmer_counts,
                      std::vector<sketch::hyperloglog> & sketches);

} // namespace seqan::hibf::sketch
