// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestPartResult, CmpHelperEQ, CmpHelperEQFailure, Assert...

#include <cstddef>    // for size_t
#include <functional> // for function
#include <ranges>     // for __fn, iota, views
#include <stdexcept>  // for runtime_error
#include <string>     // for basic_string
#include <vector>     // for vector

#include <hibf/config.hpp>                  // for config, insert_iterator
#include <hibf/sketch/compute_sketches.hpp> // for compute_sketches
#include <hibf/sketch/hyperloglog.hpp>      // for hyperloglog
#include <hibf/sketch/minhashes.hpp>        // for minhashes
#include <hibf/test/expect_range_eq.hpp>    // for expect_range_eq, EXPECT_RANGE_EQ

class compute_sketches_test : public ::testing::Test
{
public:
    std::vector<size_t> kmer_counts;
    std::vector<seqan::hibf::sketch::hyperloglog> hyperloglog_sketches;
    std::vector<seqan::hibf::sketch::minhashes> minhash_sketches;
    seqan::hibf::config config;

    void SetUp() override
    {
        config.number_of_user_bins = 3;
        config.sketch_bits = 12;
        config.input_fn = [&](size_t const ub_id, seqan::hibf::insert_iterator it)
        {
            // 0 = [0, 10000]
            // 1 = [10000, 20000]
            // 1 = [20000, 30000]
            for (size_t i = ub_id * 10000; i < (ub_id + 1) * 10000; ++i)
                it = i;
        };
    }

    void check_kmer_counts()
    {
        ASSERT_EQ(kmer_counts.size(), 3);
        EXPECT_EQ(kmer_counts[0], 18405);
        EXPECT_EQ(kmer_counts[1], 18429);
        EXPECT_EQ(kmer_counts[2], 18427);
    }

    void check_hyperloglog_sketches()
    {
        ASSERT_EQ(hyperloglog_sketches.size(), 3);
        constexpr double abs_error = 0.00001;
        EXPECT_NEAR(hyperloglog_sketches[0].estimate(), 18405.11680625227, abs_error);
        EXPECT_NEAR(hyperloglog_sketches[1].estimate(), 18429.448850770688, abs_error);
        EXPECT_NEAR(hyperloglog_sketches[2].estimate(), 18427.74236590719, abs_error);
    }

    void check_minhash_sketches()
    {
        ASSERT_EQ(minhash_sketches.size(), 3);
        size_t start{};
        size_t ub_id{};
        for (auto & minhash : minhash_sketches)
        {
            ASSERT_TRUE(minhash.is_valid()) << "ub_id: " << ub_id;
            size_t sid{};
            for (auto & sketch : minhash.table)
            {
                ASSERT_EQ(sketch.size(), seqan::hibf::sketch::minhashes::sketch_size)
                    << "ub_id: " << ub_id << "sid: " << sid;
                EXPECT_RANGE_EQ(sketch, std::views::iota(start, start + seqan::hibf::sketch::minhashes::sketch_size))
                    << "ub_id: " << ub_id << "sid: " << sid;
                ++sid;
            }
            start += 625;
            ++ub_id;
        }
    }
};

TEST_F(compute_sketches_test, hyperloglog_and_kmer_counts)
{
    seqan::hibf::sketch::compute_sketches(this->config, this->hyperloglog_sketches);

    this->check_hyperloglog_sketches();
}

TEST_F(compute_sketches_test, with_minHash)
{
    seqan::hibf::sketch::compute_sketches(this->config, this->hyperloglog_sketches, this->minhash_sketches);

    this->check_hyperloglog_sketches();
    this->check_minhash_sketches();
}

TEST_F(compute_sketches_test, too_few_hashes)
{
    this->config.number_of_user_bins = 1;
    this->config.sketch_bits = 12;
    this->config.input_fn = [&](size_t const, seqan::hibf::insert_iterator it)
    {
        for (size_t i = 0; i < 100; ++i)
            it = i;
    };

    EXPECT_THROW(
        seqan::hibf::sketch::compute_sketches(this->config, this->hyperloglog_sketches, this->minhash_sketches),
        std::runtime_error);
}
