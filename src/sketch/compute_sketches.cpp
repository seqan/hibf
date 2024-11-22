// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for __fn, sort_heap, fill, make_heap
#include <atomic>     // for atomic_flag
#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for equal_to, function
#include <limits>     // for numeric_limits
#include <span>       // for span
#include <stdexcept>  // for runtime_error
#include <string>     // for allocator, char_traits, operator+, to_string
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

/*!\brief Encapsulates handling of too few kmers to compute minHash sketches.
 * \details
 * OMP does not allow just throwing in a thread; it requires that the throw is handled in the same thread.
 * In order to properly throw the exeception in the main thread, we need to set a flag.
 * This flag then tells the other threads to stop working.
 * Afterwards, the main thread can throw the exception.
 */
struct too_few_kmers_handler
{
    std::atomic_flag flag{};
    size_t available{};

    inline bool test() noexcept
    {
        return flag.test();
    }

    inline void set(size_t const new_available) noexcept
    {
        // Sets the flag to true and returns previous value.
        // If the flag was already set, another thread encountered this block at the same time.
        // This basically acts as a mutex for setting available.
        if (!flag.test_and_set())
            this->available = new_available;
    }

    inline void check_and_throw()
    {
        if (test())
            throw std::runtime_error{"Not enough kmers (" + std::to_string(available) + ") to get "
                                     + std::to_string(minhashes::num_sketches * minhashes::sketch_size)
                                     + " hashes for all minHash sketches."};
    }
};

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

    too_few_kmers_handler too_few_kmers{};

    std::vector<uint64_t> heap{};

#pragma omp parallel for schedule(dynamic) num_threads(config.threads) private(kmers) private(heap)
    for (size_t i = 0; i < config.number_of_user_bins; ++i)
    {
        // Skip work if we already know that we have too few kmers.
        if (too_few_kmers.test())
            continue;

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

        seqan::hibf::sketch::minhashes minhash_sketch{heap};

        // In case the former heap_size wasn't sufficient to fill the minhash table
        // we need to increase the heap size and fill the heap further. We can take advantage of values
        // that already had been in the heap because they are still the smallest ones.
        while (!minhash_sketch.is_valid() && !too_few_kmers.test())
        {
            heap_size *= 2;

            if (heap_size > kmers.size())
            {
                too_few_kmers.set(kmers.size());
                break;
            }

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

    too_few_kmers.check_and_throw();
}

} // namespace seqan::hibf::sketch
