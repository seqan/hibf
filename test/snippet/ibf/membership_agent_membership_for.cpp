// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cstdint> // for uint64_t
#include <vector>  // for vector

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size
#include <hibf/misc/print.hpp>               // for print, print_t

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};
    std::vector<uint64_t> const query{126, 712, 237};

    ibf.emplace(165, seqan::hibf::bin_index{0u});
    for (auto && value : query)
    {
        ibf.emplace(value, seqan::hibf::bin_index{3u});
        ibf.emplace(value, seqan::hibf::bin_index{5u});
    }
    ibf.emplace(126, seqan::hibf::bin_index{7u});
    ibf.emplace(712, seqan::hibf::bin_index{7u});
    ibf.emplace(956, seqan::hibf::bin_index{9u});

    auto agent = ibf.membership_agent();
    // Returns all bin indices that contain at least 2 elements of the query.
    // Capture the result by reference to avoid copies.
    auto & result = agent.membership_for(query, 2u);
    seqan::hibf::print(result); // [3, 5, 7]

    // Calling `increase_bin_number_to` invalidates the agent.
    ibf.increase_bin_number_to(seqan::hibf::bin_count{60u});

    // So make sure to construct a new membership_agent.
    agent = ibf.membership_agent();
}
