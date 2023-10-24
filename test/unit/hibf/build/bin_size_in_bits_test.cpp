// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <vector> // for allocator

#include <hibf/build/bin_size_in_bits.hpp> // for bin_size_in_bits

TEST(bin_size_in_bits_test, general)
{
    EXPECT_EQ((seqan::hibf::build::bin_size_in_bits({.fpr = 0.05, .hash_count = 2, .elements = 1000})), 7903u);
}

TEST(bin_size_in_bits_test, no_elements)
{
    EXPECT_EQ((seqan::hibf::build::bin_size_in_bits({.fpr = 0.05, .hash_count = 1, .elements = 0})), 0u);
}

#ifndef NDEBUG
TEST(bin_size_in_bits_test, assert_trigger)
{
    EXPECT_DEATH((seqan::hibf::build::bin_size_in_bits({.fpr = 0.05, .hash_count = 0, .elements = 10})), "");
    EXPECT_DEATH((seqan::hibf::build::bin_size_in_bits({.fpr = 0.0, .hash_count = 2, .elements = 10})), "");
    EXPECT_DEATH((seqan::hibf::build::bin_size_in_bits({.fpr = 1.0, .hash_count = 2, .elements = 10})), "");
}
#endif
