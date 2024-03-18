// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <hibf/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an Interleaved Bloom Filter that contains 43 bins, each using 8192 bits, and 3 hash functions.
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{43u},
                                              seqan::hibf::bin_size{8192u},
                                              seqan::hibf::hash_function_count{3}};

    // Construct an Interleaved Bloom Filter that contains 43 bins, each using 256 KiBits,
    // and the default of 2 hash functions.
    seqan::hibf::interleaved_bloom_filter ibf2{seqan::hibf::bin_count{43}, seqan::hibf::bin_size{1ULL << 20}};
}
