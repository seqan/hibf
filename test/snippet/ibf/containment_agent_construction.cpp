// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <hibf/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an Interleaved Bloom Filter to be used with the containment_agent.
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{43u},
                                              seqan::hibf::bin_size{8192u},
                                              seqan::hibf::hash_function_count{3}};

    // The containment_agent can now be constructed by calling `containment_agent` on the Interleaved Bloom Filter.
    auto agent = ibf.containment_agent();

    // Calling `increase_bin_number_to` invalidates the agent.
    ibf.increase_bin_number_to(seqan::hibf::bin_count{60u});

    // So make sure to construct a new containment_agent.
    agent = ibf.containment_agent();
}
