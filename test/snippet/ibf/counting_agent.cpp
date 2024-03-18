// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cinttypes> // for uint8_t
#include <cstddef>   // for size_t
#include <ranges>    // for iota_view, operator==, __fn, iota, views
#include <vector>    // for vector

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size, hash_...
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
    auto & result = agent.bulk_count(sequence1); // Bind by `&` to avoid copies!
    seqan::hibf::print(result);                  // [20,0,0,0,10,0,0,0]

    // Search for specific values
    std::vector<size_t> const values{92, 1238, 812, 81273};
    seqan::hibf::print(agent.bulk_count(values));                      // [0,0,0,0,0,0,0,0]
    seqan::hibf::print(agent.bulk_count(std::views::iota(0u, 1024u))); // [20,0,0,0,20,0,0,10]

    // The default counters are 16 bit unsigned integer.
    // An optional template parameter can be used to specify the counter type
    auto agent2 = ibf.counting_agent<uint8_t>();
    // The returned counts are now 8 bit unsigned integers.
    seqan::hibf::print(agent2.bulk_count(sequence1)); // [20,0,0,0,10,0,0,0]
}
