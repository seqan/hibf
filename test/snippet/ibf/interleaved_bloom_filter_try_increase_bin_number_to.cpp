// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size
#include <hibf/misc/print.hpp>               // for print, print_t

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{73u}, seqan::hibf::bin_size{1024u}};
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});

    // Same bin count has no effect and returns `true`.
    bool result = ibf.try_increase_bin_number_to(seqan::hibf::bin_count{73u});
    std::cout << std::boolalpha << result << '\n'; // true
    std::cout << ibf.bin_count() << '\n';          // 73

    // Smaller bin count has no effect and returns `false`.
    result = ibf.try_increase_bin_number_to(seqan::hibf::bin_count{50u});
    std::cout << std::boolalpha << result << '\n'; // false
    std::cout << ibf.bin_count() << '\n';          // 73

    // Larger bin count and resize not required increases the bin count and returns `true`.
    result = ibf.try_increase_bin_number_to(seqan::hibf::bin_count{128u});
    std::cout << std::boolalpha << result << '\n'; // true
    std::cout << ibf.bin_count() << '\n';          // 128

    // Resize would be required, hence returns `false`.
    result = ibf.try_increase_bin_number_to(seqan::hibf::bin_count{129u});
    std::cout << std::boolalpha << result << '\n'; // false
    std::cout << ibf.bin_count() << '\n';          // 128

    // Be sure to get the agent after `try_increase_bin_number_to` as it may invalidate all agents!
    auto agent = ibf.membership_agent();

    // The content of the bins which were already present before the resize does not change
    seqan::hibf::print(agent.bulk_contains(126)); // [1,0,0,0,0,0,0,0,0,0,0,...,0]
    seqan::hibf::print(agent.bulk_contains(712)); // [0,0,0,1,0,0,0,0,0,0,0,...,0]
    seqan::hibf::print(agent.bulk_contains(237)); // [0,0,0,0,0,0,0,0,0,1,0,...,0]
}
