// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <hibf/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an Interleaved Bloom Filter to be used with the membership_agent.
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{43u},
                                              seqan::hibf::bin_size{8192u},
                                              seqan::hibf::hash_function_count{3}};

    // The membership_agent can now be constructed by calling `membership_agent` on the Interleaved Bloom Filter.
    auto agent = ibf.membership_agent();

    // Calling `increase_bin_number_to` invalidates the agent.
    ibf.increase_bin_number_to(seqan::hibf::bin_count{60u});

    // So make sure to construct a new membership_agent.
    agent = ibf.membership_agent();
}
