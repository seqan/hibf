// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __fn, sort_heap, all_of, fill, make_heap
#include <atomic>     // for atomic_flag
#include <cassert>    // for assert
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t, uint8_t
#include <functional> // for equal_to, function
#include <limits>     // for numeric_limits
#include <span>       // for span
#include <stdexcept>  // for runtime_error
#include <string>     // for operator+, allocator, char_traits, to_string, basic_string
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

        size_t heap_size{1000};
        heap.resize(heap_size);
        std::ranges::fill(heap, std::numeric_limits<uint64_t>::max()); // proper heap

        for (auto const hash : kmers)
        {
            hll_sketch.add(hash);
            seqan::hibf::sketch::minhashes::push_to_heap_if_smaller(hash, heap);
        }

        std::ranges::sort_heap(heap); // sort ascending to get the smallest numbers first

        if (heap.size() > kmers.size()) // If heap size is larger then number of available kmers to begin with
            heap.resize(kmers.size()); // such that unitializes values are cut off

        seqan::hibf::sketch::minhashes minhash_sketch{heap};

        // In case the former heap_size wasn't sufficient to fill the minhash table
        // we need to increase the heap size and fill the heap further. We can take advantage of values
        // that already had been in the heap because they are still the smallest ones.
        while (!minhash_sketch.is_valid())
        {
            heap_size *= 2;

            if (heap_size > kmers.size())
                break;

            heap.resize(heap_size, std::numeric_limits<uint64_t>::max());
            std::ranges::make_heap(heap);

            for (auto const hash : kmers)
                seqan::hibf::sketch::minhashes::push_to_heap_if_smaller(hash, heap);

            std::ranges::sort_heap(heap); // sort ascending to get the smallest numbers first

            assert(heap_size % 2u == 0u);
            std::span<uint64_t> const heap_to_consider(heap.begin() + heap_size / 2, heap.end());
            minhash_sketch.fill_incomplete_sketches(heap_to_consider);
        }

        hll_sketches[i] = std::move(hll_sketch);
        minhash_sketches[i] = std::move(minhash_sketch);
    }
}

} // namespace seqan::hibf::sketch
