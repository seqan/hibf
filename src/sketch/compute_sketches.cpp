// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __fn, sort_heap, all_of, min
#include <cassert>    // for assert
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t, uint8_t
#include <functional> // for equal_to, function
#include <limits>     // for numeric_limits
#include <span>       // for span
#include <utility>    // for move
#include <vector>     // for vector

#include <hibf/config.hpp>                  // for config, insert_iterator
#include <hibf/contrib/robin_hood.hpp>      // for hash, unordered_flat_set
#include <hibf/sketch/compute_sketches.hpp> // for compute_sketches
#include <hibf/sketch/hyperloglog.hpp>      // for hyperloglog
#include <hibf/sketch/minhashes.hpp>        // for minhashes

namespace seqan::hibf::sketch
{

void compute_sketches(config const & config, std::vector<sketch::hyperloglog> & hll_sketches)
{
    // compute hll_sketches
    hll_sketches.resize(config.number_of_user_bins, config.sketch_bits);

    assert(std::ranges::all_of(hll_sketches,
                               [bits = config.sketch_bits](hyperloglog const & sketch)
                               {
                                   return sketch.data_size() == (1ULL << bits);
                               }));

#pragma omp parallel for schedule(dynamic) num_threads(config.threads)
    for (size_t i = 0; i < config.number_of_user_bins; ++i)
    {
        config.input_fn(i, insert_iterator{hll_sketches[i]});
    }
}

// minhash_sketches data structure:
// Vector L1 : number of user bins
// Vector L2 : number_of_max_minhash_sketches (LSH ADD+OR parameter b)
// Vector L3 : minHash_sketche_size (LSH ADD+OR parameter r)
void compute_sketches(config const & config,
                      std::vector<sketch::hyperloglog> & hll_sketches,
                      std::vector<sketch::minhashes> & minhash_sketches)
{
    //inititalise
    hll_sketches.resize(config.number_of_user_bins);
    minhash_sketches.resize(config.number_of_user_bins);

    // compute hll_sketches
    robin_hood::unordered_flat_set<uint64_t> kmers{};

    std::vector<uint64_t> heap{};

#pragma omp parallel for schedule(dynamic) num_threads(config.threads) private(kmers) private(heap)
    for (size_t i = 0; i < config.number_of_user_bins; ++i)
    {
        seqan::hibf::sketch::hyperloglog hll_sketch{config.sketch_bits};

        kmers.clear();
        config.input_fn(i, insert_iterator{kmers});

        // The heap must not be larger than `kmers`; otherwise, it would keep `max()` placeholders.
        size_t heap_size = std::min<size_t>(1000u, kmers.size());
        heap.assign(heap_size, std::numeric_limits<uint64_t>::max()); // proper heap

        for (auto const hash : kmers)
        {
            hll_sketch.add(hash);
            seqan::hibf::sketch::minhashes::push_to_heap_if_smaller(hash, heap);
        }

        std::ranges::sort_heap(heap); // sort ascending to get the smallest numbers first

        seqan::hibf::sketch::minhashes minhash_sketch{heap};

        // In case the former heap_size wasn't sufficient to fill the minhash table, we need a larger heap.
        // Once sorted, it begins with the values of the previous heap, which are already in minhash_sketch.
        while (!minhash_sketch.is_valid())
        {
            // All kmers are in the heap: too few to fill the minhash table, which stays invalid.
            if (heap_size == kmers.size())
                break;

            size_t const previous_heap_size = heap_size;
            heap_size = std::min(2u * heap_size, kmers.size());

            // Start over: pushing all kmers into the previous heap would insert its values a second time.
            // Keeping the previous heap and pushing only the kmers with `wyhash(kmer) > previous_max`, its largest value,
            // saves about half of the heap operations: 15-20 % faster for up to 100'000 kmers, hardly any gain for
            // millions. It needs wyhash outside of push_to_heap_if_smaller, though.
            heap.assign(heap_size, std::numeric_limits<uint64_t>::max());

            for (auto const hash : kmers)
                seqan::hibf::sketch::minhashes::push_to_heap_if_smaller(hash, heap);

            std::ranges::sort_heap(heap); // sort ascending to get the smallest numbers first

            std::span<uint64_t> const heap_to_consider(heap.begin() + previous_heap_size, heap.end());
            minhash_sketch.fill_incomplete_sketches(heap_to_consider);
        }

        hll_sketches[i] = std::move(hll_sketch);
        minhash_sketches[i] = std::move(minhash_sketch);
    }
}

} // namespace seqan::hibf::sketch
