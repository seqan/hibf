// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size
#include <hibf/misc/print.hpp>               // for print, print_t

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});

    ibf.increase_bin_number_to(seqan::hibf::bin_count{18u});
    // Be sure to get the agent after `increase_bin_number_to` as it invalidates all agents!
    auto agent = ibf.containment_agent();

    // The content of the bins which were already present before the resize does not change
    seqan::hibf::print(agent.bulk_contains(126)); // [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    seqan::hibf::print(agent.bulk_contains(712)); // [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    seqan::hibf::print(agent.bulk_contains(237)); // [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
}
