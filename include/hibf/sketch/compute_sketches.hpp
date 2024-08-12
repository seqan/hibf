// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/config.hpp>             // for config
#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog
#include <hibf/sketch/minhashes.hpp>   // for minhash

namespace seqan::hibf::sketch
{

/*!\brief Computes the kmer_counts and sketches and stores them in the respective vectors for further use.
 * \ingroup hibf_layout
 * \param[in] config The configuration to compute the layout with.
 * \param[in,out] hll_sketches The vector that will store the sketches.
 */
void compute_sketches(config const & config, std::vector<sketch::hyperloglog> & hll_sketches);

//!\overload
void compute_sketches(config const & config,
                      std::vector<sketch::hyperloglog> & hll_sketches,
                      std::vector<sketch::minhashes> & minhash_sketches);

} // namespace seqan::hibf::sketch
