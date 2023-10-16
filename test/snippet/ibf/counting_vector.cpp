// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cinttypes> // for uint16_t

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, counting_vector, bin_count
#include <hibf/misc/print.hpp>               // for print, print_t

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(126, seqan::hibf::bin_index{3u});
    ibf.emplace(126, seqan::hibf::bin_index{9u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});

    // The counting_vector must be at least as big as the number of bins.
    seqan::hibf::counting_vector<uint16_t> counts(12, 0);

    auto agent = ibf.membership_agent();

    counts += agent.bulk_contains(712); // `counts` contains the number of occurrences of 712 in each bin.
    seqan::hibf::print(counts);         // prints [0,0,0,1,0,0,0,0,0,0,0,0]

    counts += agent.bulk_contains(237); // `counts` contains the number of occurrences of 712 and 237 in each bin.
    seqan::hibf::print(counts);         // prints [0,0,0,1,0,0,0,0,0,1,0,0]

    counts += agent.bulk_contains(126); // `counts` contains the number of occurrences of 712, 237 and 126 in each bin.
    seqan::hibf::print(counts);         // prints [1,0,0,2,0,0,0,0,0,2,0,0]

    counts += counts;           // multiple counts can also be added together
    seqan::hibf::print(counts); // prints [2,0,0,4,0,0,0,0,0,4,0,0]
}
