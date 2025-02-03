// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, TestInfo, TEST, TestPartResult, CmpHelperEQ, Cmp...

#include <algorithm>  // for __fn, copy
#include <array>      // for array
#include <cstddef>    // for size_t
#include <cstdint>    // for uint64_t, uint8_t
#include <functional> // for equal_to, function
#include <string>     // for basic_string
#include <vector>     // for vector

#include <hibf/config.hpp>                   // for insert_iterator
#include <hibf/contrib/robin_hood.hpp>       // for hash, unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_count, bin_size, hash_function_count
#include <hibf/sketch/hyperloglog.hpp>       // for hyperloglog
#include <hibf/test/expect_range_eq.hpp>     // for expect_range_eq, EXPECT_RANGE_EQ

static constexpr std::array<size_t, 10> values{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

TEST(insert_iterator_test, unordered_set)
{
    robin_hood::unordered_flat_set<uint64_t> target;
    seqan::hibf::insert_iterator it{target};
    std::ranges::copy(values, it);
    EXPECT_EQ(target.size(), 10u);
}

TEST(insert_iterator_test, sketch)
{
    seqan::hibf::sketch::hyperloglog target{5u};
    seqan::hibf::insert_iterator it{target};
    std::ranges::copy(values, it);
    EXPECT_NEAR(target.estimate(), 11.99, 0.001);
}

TEST(insert_iterator_test, ibf)
{
    seqan::hibf::interleaved_bloom_filter target{seqan::hibf::bin_count{8u},
                                                 seqan::hibf::bin_size{8u},
                                                 seqan::hibf::hash_function_count{1u}};
    for (size_t i = 0; i < 3; ++i)
    {
        seqan::hibf::insert_iterator it{target, i};
        std::ranges::copy(values, it);
    }

    auto agent = target.counting_agent<uint8_t>();
    auto & result = agent.bulk_count(values);
    std::vector<uint8_t> const expected{10, 10, 10, 0, 0, 0, 0, 0};
    EXPECT_RANGE_EQ(result, expected);
}

TEST(insert_iterator_test, function)
{
    robin_hood::unordered_flat_set<uint64_t> target;
    std::function<void(uint64_t const)> fun = [&target](size_t const value)
    {
        target.emplace(value);
        target.emplace((1u + value) * 11u);
    };
    seqan::hibf::insert_iterator it{fun};
    std::ranges::copy(values, it);
    EXPECT_EQ(target.size(), 20);
}
