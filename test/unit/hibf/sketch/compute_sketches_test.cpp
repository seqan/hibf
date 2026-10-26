// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestPartResult, CmpHelperEQ, CmpHelperEQFailure, Assert...

#include <algorithm>   // for sort
#include <array>       // for array
#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <functional>  // for function
#include <limits>      // for numeric_limits
#include <stdexcept>   // for runtime_error
#include <string>      // for basic_string
#include <string_view> // for basic_string_view
#include <utility>     // for pair
#include <vector>      // for vector

#include <hibf/config.hpp>                  // for config, insert_iterator
#include <hibf/sketch/compute_sketches.hpp> // for compute_sketches
#include <hibf/sketch/hyperloglog.hpp>      // for hyperloglog
#include <hibf/sketch/minhashes.hpp>        // for minhashes
#include <hibf/sketch/wyhash.hpp>           // for wyhash
#include <hibf/test/expect_range_eq.hpp>    // for expect_range_eq, EXPECT_RANGE_EQ

static constexpr uint64_t max_hash = std::numeric_limits<uint64_t>::max();

// The register of a hash is `hash & 15`, and minhash sketch `r` holds the 40 smallest hashes of register `r`.
// Inserts the values in [0, size), except those whose hash is in register 15 and smaller than `threshold`.
// Register 15 is then rare among the smallest hashes, and the smallest 1000 hashes do not fill sketch 15.
static void sparse_register_15(size_t const size, uint64_t const threshold, seqan::hibf::insert_iterator & it)
{
    for (uint64_t value = 0; value < size; ++value)
        if (uint64_t const hash = seqan::hibf::sketch::wyhash(value); (hash & 15u) != 15u || hash >= threshold)
            it = value;
}

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
        ASSERT_EQ(minhash_sketches.size(), config.number_of_user_bins);
        size_t ub_id{};
        for (auto & minhash : minhash_sketches)
        {
            ASSERT_TRUE(minhash.is_valid()) << "ub_id: " << ub_id;
            std::vector<uint64_t> hashes{};
            seqan::hibf::insert_iterator::function_t push_hash = [&hashes](uint64_t const value)
            {
                hashes.push_back(seqan::hibf::sketch::wyhash(value));
            };
            config.input_fn(ub_id, seqan::hibf::insert_iterator{push_hash});
            std::ranges::sort(hashes);
            seqan::hibf::sketch::minhashes const expected{hashes};
            size_t sid{};
            for (auto & sketch : minhash.table)
            {
                ASSERT_EQ(sketch.hashes.size(), seqan::hibf::sketch::minhashes::sketch_size)
                    << "ub_id: " << ub_id << "sid: " << sid;
                EXPECT_RANGE_EQ(sketch.hashes, expected.table[sid].hashes) << "ub_id: " << ub_id << "sid: " << sid;
                ++sid;
            }
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
        for (size_t i = 0; i < 10; ++i)
            it = i;
    };

    seqan::hibf::sketch::compute_sketches(this->config, this->hyperloglog_sketches, this->minhash_sketches);

    ASSERT_EQ(this->minhash_sketches.size(), 1); // only one user bin

    for (auto const & sketch : this->minhash_sketches[0].table)
    {
        ASSERT_FALSE(sketch.is_valid()) << "It by chance is valid";
    }
}

TEST_F(compute_sketches_test, minHash_refill)
{
    // compute_sketches first keeps the smallest 1000 hashes of a user bin. If they do not fill all 16 sketches, it
    // refills with a larger heap (2000, 4000, ..., at most the number of kmers). Here, both user bins need a refill:
    // * User bin 0: 3968 kmers. Register 15 has 33 of the smallest 1000 and 95 of the smallest 2000 hashes.
    //   The old refill pushed all kmers into the previous heap. The heap of 2000 then held each of the smallest 1000
    //   hashes twice, and sketch 15 got duplicates.
    // * User bin 1: 1463 kmers. Register 15 has 26 of the smallest 1000 hashes and 56 in total. The heap has to grow
    //   to 1463. The old code threw instead, because 2000 > 1463.
    // check_minhash_sketches compares the result to the sketches built from all hashes at once.
    // Each user bin: {number of values, threshold}.
    static constexpr std::array<std::pair<size_t, uint64_t>, 2> bins{
        {{4000u, max_hash / 8u}, {1500u, max_hash / 5u * 2u}}};
    this->config.number_of_user_bins = bins.size();
    this->config.input_fn = [](size_t const ub_id, seqan::hibf::insert_iterator it)
    {
        sparse_register_15(bins[ub_id].first, bins[ub_id].second, it);
    };

    seqan::hibf::sketch::compute_sketches(this->config, this->hyperloglog_sketches, this->minhash_sketches);

    this->check_minhash_sketches();
}

TEST_F(compute_sketches_test, too_few_hashes_in_one_register)
{
    // 960 kmers, but only 31 hashes in register 15, so sketch 15 cannot be filled and stays incomplete.
    // With fewer than 1000 kmers, the old code kept max() placeholders in the heap. These filled sketch 15 with
    // max() >> 4.
    this->config.number_of_user_bins = 1;
    this->config.input_fn = [](size_t const, seqan::hibf::insert_iterator it)
    {
        sparse_register_15(990u, max_hash / 2u, it);
    };

    seqan::hibf::sketch::compute_sketches(this->config, this->hyperloglog_sketches, this->minhash_sketches);

    ASSERT_EQ(this->minhash_sketches.size(), 1u);
    auto const & table = this->minhash_sketches[0].table;
    ASSERT_EQ(table.size(), seqan::hibf::sketch::minhashes::num_sketches);
    EXPECT_FALSE(this->minhash_sketches[0].is_valid());
    EXPECT_EQ(table[15].hashes.size(), 31u); // no max() >> 4 placeholders
}
