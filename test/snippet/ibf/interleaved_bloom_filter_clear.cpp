// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges> // for iota_view, operator==, _Iota, iota, views
#include <vector> // for vector

#include <hibf/interleaved_bloom_filter.hpp> // for bin_index, interleaved_bloom_filter, bin_count, bin_size, hash_...
#include <hibf/misc/print.hpp>               // for print, print_t

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{8u},
                                              seqan::hibf::bin_size{8192u},
                                              seqan::hibf::hash_function_count{2u}};

    auto sequence1 = std::views::iota(0u, 20u);
    auto sequence2 = std::views::iota(10u, 30u);
    auto sequence3 = std::views::iota(25u, 35u);

    // Insert all values of sequence1 into bin 0
    for (auto && value : sequence1)
        ibf.emplace(value, seqan::hibf::bin_index{0u});

    // Insert all values of sequence2 into bin 4
    for (auto && value : sequence2)
        ibf.emplace(value, seqan::hibf::bin_index{4u});

    // Insert all values of sequence3 into bin 7
    for (auto && value : sequence3)
        ibf.emplace(value, seqan::hibf::bin_index{7u});

    auto agent = ibf.counting_agent();

    // Count all values of sequence1 for all bins
    seqan::hibf::print(agent.bulk_count(sequence1)); // [20,0,0,0,10,0,0,0]

    // Clear bin 0
    ibf.clear(seqan::hibf::bin_index{0u});

    // After clearing, no values are found in bin 0
    seqan::hibf::print(agent.bulk_count(sequence1)); // [0,0,0,0,10,0,0,0]

    // Search for specific values
    seqan::hibf::print(agent.bulk_count(std::views::iota(0u, 1024u))); // [0,0,0,0,20,0,0,10]

    // Clear bin 4 and 7
    ibf.clear(std::vector{seqan::hibf::bin_index{4u}, seqan::hibf::bin_index{7u}});

    // After clearing, nothing is found
    seqan::hibf::print(agent.bulk_count(std::views::iota(0u, 1024u))); // [0,0,0,0,0,0,0,0]
}
