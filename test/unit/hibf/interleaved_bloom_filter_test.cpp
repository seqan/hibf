// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HIBF_HAS_AVX512
#    define HIBF_HAS_AVX512 0
#endif

#include <gtest/gtest.h> // for Message, Test, AssertionResult, TestPartResult, TestInfo, EXPEC...

#include <algorithm>   // for __fn, for_each
#include <array>       // for array
#include <cinttypes>   // for uint64_t
#include <compare>     // for operator<, strong_ordering
#include <cstddef>     // for size_t
#include <functional>  // for function
#include <ranges>      // for iota_view, __fn, iota, views, operator==
#include <stdexcept>   // for logic_error, invalid_argument
#include <type_traits> // for is_copy_assignable_v, is_copy_constructible_v, is_default_const...
#include <utility>     // for move
#include <vector>      // for vector

#include <hibf/config.hpp>                   // for insert_iterator, config
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size, hash_...
#include <hibf/misc/bit_vector.hpp>          // for bit_vector
#include <hibf/misc/counting_vector.hpp>     // for counting_vector
#include <hibf/test/cereal.hpp>              // for test_serialisation
#include <hibf/test/expect_range_eq.hpp>     // for expect_range_eq, EXPECT_RANGE_EQ

TEST(ibf_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<seqan::hibf::interleaved_bloom_filter>);
    EXPECT_TRUE(std::is_copy_constructible_v<seqan::hibf::interleaved_bloom_filter>);
    EXPECT_TRUE(std::is_move_constructible_v<seqan::hibf::interleaved_bloom_filter>);
    EXPECT_TRUE(std::is_copy_assignable_v<seqan::hibf::interleaved_bloom_filter>);
    EXPECT_TRUE(std::is_move_assignable_v<seqan::hibf::interleaved_bloom_filter>);
    EXPECT_TRUE(std::is_destructible_v<seqan::hibf::interleaved_bloom_filter>);

    // num hash functions defaults to two
    seqan::hibf::interleaved_bloom_filter ibf1{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{1024u}};
    seqan::hibf::interleaved_bloom_filter ibf2{seqan::hibf::bin_count{64u},
                                               seqan::hibf::bin_size{1024u},
                                               seqan::hibf::hash_function_count{2u}};
    EXPECT_TRUE(ibf1 == ibf2);

    // bin_size parameter is too small
    EXPECT_THROW((seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{0u}}),
                 std::logic_error);
    // not enough bins
    EXPECT_THROW((seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{0u}, seqan::hibf::bin_size{32u}}),
                 std::logic_error);
    // not enough hash functions
    EXPECT_THROW((seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{64u},
                                                        seqan::hibf::bin_size{32u},
                                                        seqan::hibf::hash_function_count{0u}}),
                 std::logic_error);
    // too many hash functions
    EXPECT_THROW((seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{64u},
                                                        seqan::hibf::bin_size{32u},
                                                        seqan::hibf::hash_function_count{6u}}),
                 std::logic_error);
}

TEST(ibf_test, construction_from_config)
{
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {0u, 2u, 3u, 4u, 5u}};
    size_t const number_of_ub{hashes.size()};

    seqan::hibf::config ibf_config{.input_fn =
                                       [&](size_t const num, seqan::hibf::insert_iterator it)
                                   {
                                       for (auto const hash : hashes[num])
                                           it = hash;
                                   },
                                   .number_of_user_bins = number_of_ub};

    seqan::hibf::interleaved_bloom_filter ibf{ibf_config};

    auto agent = ibf.membership_agent();

    std::vector<size_t> query{1, 2, 3, 4, 5};

    // value 2 is in both user bins
    std::vector<bool> expected_v2(number_of_ub);
    expected_v2[0] = 1;
    expected_v2[1] = 1;
    // value 8 is only in user bin 0
    std::vector<bool> expected_v8(number_of_ub);
    expected_v8[0] = 1;
    // value 0 is only in user bin 1
    std::vector<bool> expected_v0(number_of_ub);
    expected_v0[1] = 1;

    EXPECT_RANGE_EQ(agent.bulk_contains(2), expected_v2);
    EXPECT_RANGE_EQ(agent.bulk_contains(8), expected_v8);
    EXPECT_RANGE_EQ(agent.bulk_contains(0), expected_v0);
}

TEST(ibf_test, construction_from_config_with_max_bin_elements)
{
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {0u, 2u, 3u, 4u, 5u}};
    size_t const number_of_ub{hashes.size()};

    seqan::hibf::config ibf_config{.input_fn =
                                       [&](size_t const num, seqan::hibf::insert_iterator it)
                                   {
                                       for (auto const hash : hashes[num])
                                           it = hash;
                                   },
                                   .number_of_user_bins = number_of_ub};

    seqan::hibf::interleaved_bloom_filter only_config{ibf_config};
    seqan::hibf::interleaved_bloom_filter default_num_elements{ibf_config, 0u};
    seqan::hibf::interleaved_bloom_filter appropriate_num_elements{ibf_config, 10u};
    seqan::hibf::interleaved_bloom_filter larger_num_elements{ibf_config, 20u};

    EXPECT_EQ(only_config, default_num_elements);
    EXPECT_EQ(only_config, appropriate_num_elements);
    EXPECT_NE(only_config, larger_num_elements);

    EXPECT_EQ(default_num_elements, appropriate_num_elements);
    EXPECT_NE(default_num_elements, larger_num_elements);

    EXPECT_NE(appropriate_num_elements, larger_num_elements);
}

TEST(ibf_test, member_getter)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{1024u}};
    EXPECT_EQ(ibf.bin_count(), 64u);
    EXPECT_EQ(ibf.bin_size(), 1024u);
    EXPECT_EQ(ibf.bit_size(), 65'536ull);
    EXPECT_EQ(ibf.hash_function_count(), 2u);

    ibf = seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{73u},
                                                seqan::hibf::bin_size{1019u},
                                                seqan::hibf::hash_function_count{3u}};
    EXPECT_EQ(ibf.bin_count(), 73u);
    EXPECT_EQ(ibf.bin_size(), 1019u);
    EXPECT_EQ(ibf.bit_size(), 130'432ull);
    EXPECT_EQ(ibf.hash_function_count(), 3u);
}

TEST(ibf_test, bulk_contains)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{1024u}};
    std::vector<bool> expected(64); // empty bitvector is expected since we did not insert anything
    auto agent = ibf.membership_agent();

    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }

    // Test iterator interface
    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        size_t i = 0;
        for (auto it = res.begin(); it < res.end(); ++it, ++i)
        {
            EXPECT_EQ(*it, expected[i]);
        }
        EXPECT_EQ(i, expected.size());
    }

    // Test operator[] interface
    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_EQ(expected.size(), res.size());
        for (size_t i = 0; i < res.size(); ++i)
        {
            EXPECT_EQ(res[i], expected[i]);
        }
    }
}

TEST(ibf_test, emplace)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Test for correctness
    auto agent = ibf.membership_agent();
    std::vector<bool> expected(64, 1);          // every hash value should be set for every bin
    for (size_t hash : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }
}

TEST(ibf_test, emplace_exists)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{128u},
                                              seqan::hibf::bin_size{512},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Test for correctness
    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ASSERT_TRUE(ibf.emplace_exists(hash, seqan::hibf::bin_index{bin_idx}));

    for (size_t bin_idx : std::views::iota(64, 128))
        ASSERT_FALSE(ibf.emplace_exists(0u, seqan::hibf::bin_index{bin_idx}));
}

TEST(ibf_test, clear)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Clear a bin
    ibf.clear(seqan::hibf::bin_index{17u});

    // 3. Test for correctness
    auto agent = ibf.membership_agent();
    std::vector<bool> expected(64, 1); // every hash value should be set for every bin...
    expected[17] = 0;                  // ...except bin 17
    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }
}

TEST(ibf_test, clear_range)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        for (size_t hash : std::views::iota(0, 64))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Clear a range of bins
    std::vector<seqan::hibf::bin_index> bin_range{seqan::hibf::bin_index{8u},
                                                  seqan::hibf::bin_index{17u},
                                                  seqan::hibf::bin_index{45u}};
    ibf.clear(bin_range);

    // 3. Test for correctness
    auto agent = ibf.membership_agent();
    std::vector<bool> expected(64, 1); // every hash value should be set for every bin...
    expected[8] = 0;                   // ...except bin 8
    expected[17] = 0;                  // ...except bin 17
    expected[45] = 0;                  // ...except bin 45
    for (size_t hash : std::views::iota(0, 64))
    {
        auto & res = agent.bulk_contains(hash);
        EXPECT_RANGE_EQ(res, expected);
    }
}

TEST(ibf_test, counting)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{128u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 128))
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Test for correctness
    seqan::hibf::counting_vector<size_t> counting(128, 0);
    auto agent = ibf.membership_agent();
    for (size_t hash : std::views::iota(0, 128)) // test correct resize for each bin individually
    {
        counting += agent.bulk_contains(hash);
    }
    std::vector<size_t> expected(128, 128);
    EXPECT_RANGE_EQ(counting, expected);

    // Counting vectors can be added together
    std::vector<size_t> expected2(128, 256);
    counting += counting;
    EXPECT_RANGE_EQ(counting, expected2);

    // minus bit_vector
    for (size_t hash : std::views::iota(0, 128)) // test correct resize for each bin individually
    {
        counting -= agent.bulk_contains(hash);
    }
    EXPECT_RANGE_EQ(counting, std::vector<size_t>(128, 128));

    // minus other counting vector
    counting -= seqan::hibf::counting_vector<size_t>(128, 128 - 42);
    EXPECT_RANGE_EQ(counting, std::vector<size_t>(128, 42));
}

TEST(ibf_test, counting_agent)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{128u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 128))
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Test for correctness
    auto agent = ibf.counting_agent();
    auto agent2 = ibf.template counting_agent<size_t>();

    std::vector<size_t> expected(128, 128);
    EXPECT_RANGE_EQ(agent.bulk_count(std::views::iota(0u, 128u)), expected);
    EXPECT_RANGE_EQ(agent2.bulk_count(std::views::iota(0u, 128u)), expected);
}

// Check special case where there is only one `1` in the bitvector.
TEST(ibf_test, counting_no_ub)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{128u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::array<size_t, 2>{63, 127})
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Test for correctness
    seqan::hibf::counting_vector<size_t> counting(128, 0);
    auto agent = ibf.membership_agent();
    for (size_t hash : std::views::iota(0, 128)) // test correct resize for each bin individually
    {
        counting += agent.bulk_contains(hash);
    }
    std::vector<size_t> expected(128, 0);
    expected[63] = 128;
    expected[127] = 128;
    EXPECT_RANGE_EQ(counting, expected);

    // Counting vectors can be added together
    std::vector<size_t> expected2(128, 0);
    expected2[63] = 256;
    expected2[127] = 256;
    counting += counting;
    EXPECT_RANGE_EQ(counting, expected2);
}

// Check special case where there is only one `1` in the bitvector.
// Also checks that counting with `b` bins and `b % 64 != 0` works.
TEST(ibf_test, counting_agent_no_ub)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{127u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::array<size_t, 2>{63, 126})
        for (size_t hash : std::views::iota(0, 128))
            ibf.emplace(hash, seqan::hibf::bin_index{bin_idx});

    // 2. Test for correctness
    auto agent = ibf.counting_agent();
    auto agent2 = ibf.template counting_agent<size_t>();

    std::vector<size_t> expected(127, 0);
    expected[63] = 128;
    expected[126] = 128;
    EXPECT_RANGE_EQ(agent.bulk_count(std::views::iota(0u, 128u)), expected);
    EXPECT_RANGE_EQ(agent2.bulk_count(std::views::iota(0u, 128u)), expected);
}

TEST(ibf_test, increase_bin_number_to)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{73u}, seqan::hibf::bin_size{1024u}};
    size_t const original_bitsize{ibf.bit_size()};

    // 1. Throw if trying to reduce number of bins.
    EXPECT_THROW(ibf.increase_bin_number_to(seqan::hibf::bin_count{62u}), std::invalid_argument);
    EXPECT_EQ(ibf.bin_count(), 73u);
    EXPECT_EQ(ibf.bit_size(), original_bitsize);

    // 2. No change in bin_words implies no change in size.
    ibf.increase_bin_number_to({seqan::hibf::bin_count{127u}});
    EXPECT_EQ(ibf.bit_size(), original_bitsize);
    EXPECT_EQ(ibf.bin_count(), 127u);

    // 3. If resizing takes place, the inserted values must still be valid.
    auto hashes = std::views::iota(0, 64);
    for (size_t current_bin : std::views::iota(0, 64)) // test correct resize for each bin individually
    {
        seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{8u}};
        EXPECT_EQ(ibf.bit_size(), 512u);
        std::ranges::for_each(hashes,
                              [&ibf, &current_bin](auto const h)
                              {
                                  ibf.emplace(h, seqan::hibf::bin_index{current_bin});
                              });

        ibf.increase_bin_number_to(seqan::hibf::bin_count{73u});

        EXPECT_EQ(ibf.bin_count(), 73u);
        EXPECT_EQ(ibf.bit_size(), 1024u);

        std::vector<bool> expected(73, 0);
        expected[current_bin] = 1; // none of the bins except current_bin stores the hash values.
        auto agent = ibf.membership_agent();
        for (size_t const h : hashes)
        {
            auto & res = agent.bulk_contains(h);
            EXPECT_RANGE_EQ(res, expected);
        }
    }
}

TEST(ibf_test, copy_agents)
{
    // 1. Construct and emplace
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    for (size_t bin_idx : std::views::iota(0, 64))
        ibf.emplace(0u, seqan::hibf::bin_index{bin_idx});

    {
        auto membership_agent1 = ibf.membership_agent();
        auto membership_agent2 = membership_agent1;

        auto & result1 = membership_agent1.bulk_contains(0u);
        auto & result2 = membership_agent2.bulk_contains(1u);
        ASSERT_EQ(result1.size(), result2.size());

        EXPECT_TRUE(result1.all());
        EXPECT_TRUE(result2.none());
    }
    {
        auto counting_agent1 = ibf.counting_agent();
        auto counting_agent2 = counting_agent1;
        auto & result1 = counting_agent1.bulk_count(std::vector<uint64_t>{0u});
        auto & result2 = counting_agent2.bulk_count(std::vector<uint64_t>{1u});
        ASSERT_EQ(result1.size(), result2.size());
        EXPECT_RANGE_EQ(result1, std::vector<size_t>(64u, 1u));
        EXPECT_RANGE_EQ(result2, std::vector<size_t>(64u, 0u));
    }
}

TEST(ibf_test, serialisation)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{128u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{3u}};
    seqan::hibf::test::test_serialisation(std::move(ibf));
}
