// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, CmpHelperEQ, CmpHelperEQFailure, TestPartResult

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t
#include <random>  // for uniform_int_distribution, mt19937_64
#include <string>  // for basic_string
#include <vector>  // for vector

#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog

TEST(estimate_kmer_counts_test, small_example)
{
    seqan::hibf::sketch::hyperloglog sketch(12);

    std::uniform_int_distribution<uint64_t> distribution{};
    std::mt19937_64 engine{0u};

    for (size_t i = 0; i < 1500; ++i)
        sketch.add(distribution(engine));

    std::vector<seqan::hibf::sketch::hyperloglog> sketches{sketch, sketch};
    std::vector<size_t> kmer_counts;

    seqan::hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);

    ASSERT_EQ(kmer_counts.size(), 2);
    EXPECT_EQ(kmer_counts[0], 1498);
    EXPECT_EQ(kmer_counts[1], 1498);
}
