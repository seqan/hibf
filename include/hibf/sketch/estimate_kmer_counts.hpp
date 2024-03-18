// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

namespace seqan::hibf::sketch
{

/*!\brief Estimates k-mer counts via sketches.
 * \ingroup hibf_sketch
 */
inline void estimate_kmer_counts(std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                 std::vector<size_t> & kmer_counts)
{
    kmer_counts.resize(sketches.size());

    for (size_t i = 0; i < sketches.size(); ++i)
        kmer_counts[i] = sketches[i].estimate();
}

} // namespace seqan::hibf::sketch
