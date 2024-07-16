// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>
#include <cassert>

#include <hibf/sketch/minhashes.hpp>

namespace seqan::hibf::sketch
{

minhashes::minhashes(std::vector<uint64_t> const & smallest_values)
{
    assert(std::ranges::is_sorted(smallest_values));

    table.resize(num_sketches);

    for (auto & elem : table)
        elem.reserve(sketch_size);

    for (uint64_t const hash : smallest_values)
    {
        auto & hash_table = table[hash & register_id_mask];
        if (hash_table.size() < sketch_size)
            hash_table.push_back(hash >> 4);
    }
}

bool minhashes::is_valid() const
{
    return table.size() == num_sketches
        && std::ranges::all_of(table,
                               [](auto const & minHash_sketch)
                               {
                                   return minHash_sketch.size() == sketch_size;
                               });
}

void minhashes::fill_incomplete_sketches(std::span<uint64_t> const & more_smallest_values)
{
    assert(std::ranges::is_sorted(more_smallest_values));

    for (uint64_t const hash : more_smallest_values)
    {
        auto & hash_table = table[hash & register_id_mask];
        assert(std::ranges::find(hash_table, hash) == hash_table.end()); // hashes should be unique
        if (hash_table.size() < sketch_size)
            hash_table.push_back(hash >> 4);
    }
}

void minhashes::push_to_heap_if_smaller(uint64_t const value, std::vector<uint64_t> & heap)
{
    // Do nothing if value is bigger than the current biggest element in the (max) heap.
    if (value >= heap[0])
        return;

    // we do not need a guard (hash table) to check for duplications because `kmers` is already a set
    std::ranges::pop_heap(heap);  // max elements move to end of vector
    heap.back() = value;          // replace last elements instead of physically popping and pushing
    std::ranges::push_heap(heap); // last elements is rearranged in the heap to be pushed
}

} // namespace seqan::hibf::sketch
