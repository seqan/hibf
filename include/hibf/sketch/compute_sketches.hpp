// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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

void compute_sketches(config const & config,
                      std::vector<std::vector<size_t>> & kmer_counts,
                      std::vector<std::vector<sketch::hyperloglog>> & sketches,
                      size_t const number_of_partitions);

} // namespace seqan::hibf::sketch
