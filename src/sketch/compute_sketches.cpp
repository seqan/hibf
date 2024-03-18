// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <vector>     // for vector

#include <hibf/config.hpp>                      // for config, insert_iterator
#include <hibf/contrib/robin_hood.hpp>          // for unordered_flat_set
#include <hibf/sketch/compute_sketches.hpp>     // for compute_sketches
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog

namespace seqan::hibf::sketch
{

void compute_sketches(config const & config,
                      std::vector<size_t> & kmer_counts,
                      std::vector<sketch::hyperloglog> & sketches)
{
    // compute sketches
    sketches.resize(config.number_of_user_bins);
    kmer_counts.resize(config.number_of_user_bins);

    robin_hood::unordered_flat_set<uint64_t> kmers;
#pragma omp parallel for schedule(dynamic) num_threads(config.threads) private(kmers)
    for (size_t i = 0; i < config.number_of_user_bins; ++i)
    {
        seqan::hibf::sketch::hyperloglog sketch(config.sketch_bits);

        kmers.clear();
        config.input_fn(i, insert_iterator{kmers});

        for (auto k_hash : kmers)
            sketch.add(k_hash);

        // #pragma omp critical
        sketches[i] = sketch;
    }

    sketch::estimate_kmer_counts(sketches, kmer_counts);
}

} // namespace seqan::hibf::sketch
