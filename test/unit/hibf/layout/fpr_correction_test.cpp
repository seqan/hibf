// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, TestPartResult, EXPECT_NEAR, TestInfo

#include <cstddef> // for size_t
#include <vector>  // for vector, allocator

#include <hibf/layout/compute_fpr_correction.hpp> // for compute_fpr_correction

TEST(fp_correction_test, one_bin)
{

    auto fp_correction = seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2u, .t_max = 8u});

    std::vector<size_t> const values{9123, 123, 12, 87123, 8123, 4660};

    // Splitting into 1 bin, i.e. not splitting, should not change the bin size.
    for (size_t const value : values)
        EXPECT_EQ(value, value * fp_correction[1]);
}

TEST(fp_correction_test, example_split)
{
    auto fp_correction = seqan::hibf::layout::compute_fpr_correction({.fpr = 0.01, .hash_count = 5u, .t_max = 256u});

    double const abs_error = 0.00001;
    EXPECT_NEAR(fp_correction[1], 1.0, abs_error);
    EXPECT_NEAR(fp_correction[2], 1.192316, abs_error);
    EXPECT_NEAR(fp_correction[4], 1.412390, abs_error);
    EXPECT_NEAR(fp_correction[8], 1.664459, abs_error);
    EXPECT_NEAR(fp_correction[16], 1.953384, abs_error);
    EXPECT_NEAR(fp_correction[32], 2.284738, abs_error);
    EXPECT_NEAR(fp_correction[64], 2.664909, abs_error);
    EXPECT_NEAR(fp_correction[128], 3.101225, abs_error);
    EXPECT_NEAR(fp_correction[256], 3.602093, abs_error);

    ASSERT_EQ(fp_correction.size(), 257u);
    for (size_t i{1u}; i < 256u; ++i)
        ASSERT_LE(fp_correction[i], fp_correction[i + 1u]);
}
