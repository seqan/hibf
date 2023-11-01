// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, TestPartResult, TestInfo, EXPECT_EQ

#include <cinttypes>  // for uint32_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <ranges>     // for _Iota, iota, views
#include <sstream>    // for basic_stringstream, char_traits, stringstream
#include <utility>    // for move
#include <vector>     // for vector, allocator

#include <hibf/config.hpp>                                // for insert_iterator, config
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/interleaved_bloom_filter.hpp>              // for counting_vector
#include <hibf/layout/layout.hpp>                         // for layout
#include <hibf/test/cereal.hpp>                           // for test_serialisation
#include <hibf/test/expect_range_eq.hpp>                  // for expect_range_eq, EXPECT_RANGE_EQ

TEST(hibf_test, small_example_with_direct_hashes)
{
    // range of range of sequences
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    seqan::hibf::config config{.input_fn =
                                   [&](size_t const num, seqan::hibf::insert_iterator it)
                               {
                                   for (auto const hash : hashes[num])
                                       it = hash;
                               },
                               .number_of_user_bins = 2};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    auto agent = hibf.membership_agent();

    std::vector<size_t> query{1u, 2u, 3u, 4u, 5u};

    auto & result = agent.membership_for(query, 2u);
    agent.sort_results();
    EXPECT_RANGE_EQ(result, (std::vector<size_t>{0u, 1u}));
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);
}

TEST(hibf_test, simple_user_bins)
{
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const, seqan::hibf::insert_iterator it)
                               {
                                   it = 5u;
                               },
                               .number_of_user_bins = 73};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);
}

TEST(hibf_test, build_from_layout)
{
    // range of range of sequences
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    auto input_fn = [&](size_t const num, seqan::hibf::insert_iterator it)
    {
        for (auto const hash : hashes[num])
            it = hash;
    };

    std::stringstream stream{"@HIBF_CONFIG\n"
                             "@{\n"
                             "@    \"hibf_config\": {\n"
                             "@        \"version\": 1,\n"
                             "@        \"number_of_user_bins\": 2,\n"
                             "@        \"number_of_hash_functions\": 2,\n"
                             "@        \"maximum_fpr\": 0.05,\n"
                             "@        \"relaxed_fpr\": 0.3,\n"
                             "@        \"threads\": 1,\n"
                             "@        \"sketch_bits\": 12,\n"
                             "@        \"tmax\": 64,\n"
                             "@        \"alpha\": 1.2,\n"
                             "@        \"max_rearrangement_ratio\": 0.5,\n"
                             "@        \"disable_estimate_union\": false,\n"
                             "@        \"disable_rearrangement\": true\n"
                             "@    }\n"
                             "@}\n"
                             "@HIBF_CONFIG_END\n"
                             "#TOP_LEVEL_IBF fullest_technical_bin_idx:0\n"
                             "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n"
                             "1\t0\t34\n"
                             "0\t34\t30\n"};

    seqan::hibf::config configuration{};
    configuration.read_from(stream);
    configuration.input_fn = input_fn;

    seqan::hibf::layout::layout layout{};
    layout.read_from(stream);

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{configuration, layout};
    auto agent = hibf.membership_agent();

    std::vector<size_t> query{1u, 2u, 3u, 4u, 5u};

    auto & result = agent.membership_for(query, 2u);
    agent.sort_results();
    EXPECT_RANGE_EQ(result, (std::vector<size_t>{0u, 1u}));
    EXPECT_EQ(configuration.number_of_user_bins, hibf.number_of_user_bins);
}

TEST(hibf_test, three_level_hibf)
{
    // To ensure that there are 3 levels, we generate 4097 user bins, equal in size, with little overlap/similarity.
    // Disable rearrangement as we are not testing the quality of the layout but just that an index with three level
    // works as expected.
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const ub_id, seqan::hibf::insert_iterator it)
                               {
                                   // start every 16 positions and take 20 values
                                   // -> neighbouring user bins have an overlap of 4 hashes
                                   size_t const start = ub_id * 16;
                                   for (size_t i = start; i < start + 20; ++i)
                                       it = i;
                               },
                               .number_of_user_bins = 4097,
                               .maximum_fpr = 0.001,
                               .threads = 4,
                               .tmax = 64,
                               .disable_estimate_union = true,
                               .disable_rearrangement = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    auto agent = hibf.membership_agent();

    // test a query within a single user bin and one in the overlap
    for (size_t ub_id = 0u; ub_id < 4096u; ub_id += 100u)
    {
        size_t const start = ub_id * 16u;
        auto overlap_query = std::views::iota(start + 16u, start + 20u);
        auto unique_query = std::views::iota(start + 5u, start + 9u);

        auto & overlap_result = agent.membership_for(overlap_query, 4u); // one overlapping user bin
        agent.sort_results();
        EXPECT_EQ(overlap_result.size(), 2u);
        EXPECT_EQ(overlap_result[0], ub_id);
        EXPECT_EQ(overlap_result[1], ub_id + 1u);

        auto & unique_result = agent.membership_for(unique_query, 4u); // no overlapping user bin
        EXPECT_EQ(unique_result.size(), 1u);
        EXPECT_EQ(unique_result[0], ub_id);
    }
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);
}

// std::pow(number, 2) does conversion to double and back to size_t.
// If used more often, consider porting pow with integer overloads from
// https://github.com/seqan/seqan3/blob/8b7d02cb8695369b7baeb4b3042ae7c864b67b8c/include/seqan3/utility/math.hpp
inline size_t squared(size_t const number)
{
    return number * number;
}

TEST(hibf_test, unevenly_sized_and_unique_user_bins)
{
    // Simulate 500 user bins of exponentially increasing size
    // No overlap between the user bins happens.
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const ub_id, seqan::hibf::insert_iterator it)
                               {
                                   size_t const start = squared(ub_id + 1u);
                                   size_t const end = squared(ub_id + 2u) - 1u;
                                   for (size_t i = start; i < end; ++i)
                                       it = i;
                               },
                               .number_of_user_bins = 500,
                               .maximum_fpr = 0.001,
                               .threads = 4,
                               .tmax = 64,
                               .disable_estimate_union = true,
                               .disable_rearrangement = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    auto agent = hibf.membership_agent();

    for (size_t ub_id = 5u /* don't test the tiny ones */; ub_id < 499u; ub_id += 20u)
    {
        size_t const start = squared(ub_id + 1u);
        auto unique_query = std::views::iota(start + 5u, start + 9u);

        auto & unique_result = agent.membership_for(unique_query, 4u);
        EXPECT_EQ(unique_result.size(), 1u);
        EXPECT_EQ(unique_result[0], ub_id);
    }
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);
}

TEST(hibf_test, evenly_sized_and_highly_similar_user_bins)
{
    // Simulate 1000 user bins of roughly the same size with a lot of overlap/similarity.
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const ub_id, seqan::hibf::insert_iterator it)
                               {
                                   // Each user bin has an overlap of 30 with its predecessor,
                                   // 25 with its second predecessor, etc.
                                   // There are in total 6 other user bins that share at least 5 values.
                                   size_t const start = ub_id * 5u;
                                   for (size_t i = start; i < start + 35u; ++i)
                                       it = i;
                               },
                               .number_of_user_bins = 1000,
                               .maximum_fpr = 0.001,
                               .threads = 4,
                               .tmax = 64,
                               .disable_estimate_union = true,
                               .disable_rearrangement = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    auto agent = hibf.membership_agent();

    for (size_t ub_id = 6u; ub_id < 994u; ub_id += 20u)
    {
        size_t const start = ub_id * 5;
        auto similar_query = std::views::iota(start, start + 5u);

        auto & similar_result = agent.membership_for(similar_query, 5u); // t = 5 results in 6 similar user bins
        agent.sort_results();
        EXPECT_RANGE_EQ(similar_result, (std::views::iota(ub_id - 6u, ub_id + 1u)));
    }
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);
}

TEST(hibf_test, counting_agent_same_bins)
{
    // Simulate 500 bins, with the same 500 elements each
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const, seqan::hibf::insert_iterator it)
                               {
                                   for (size_t i = 0u; i < 500u; ++i)
                                       it = i;
                               },
                               .number_of_user_bins = 500,
                               .maximum_fpr = 0.001,
                               .threads = 4,
                               .tmax = 64,
                               .disable_estimate_union = true,
                               .disable_rearrangement = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    auto agent = hibf.template counting_agent<uint32_t>();

    auto query = std::views::iota(0u, 500u);
    auto & result = agent.bulk_count(query, 1u);

    for (size_t ub_id = 0u; ub_id < config.number_of_user_bins; ++ub_id)
    {
        ASSERT_GE(result[ub_id], 500u) << "ub_id: " << ub_id;
        // Split bins might have FP
        ASSERT_LE(result[ub_id], 502u) << "ub_id: " << ub_id;
    }
}

TEST(hibf_test, counting_agent_different_bins)
{
    // Simulate 500 user bins of exponentially increasing size
    // No overlap between the user bins happens.
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const ub_id, seqan::hibf::insert_iterator it)
                               {
                                   size_t const start = squared(ub_id + 1u);
                                   size_t const end = squared(ub_id + 2u) - 1u;
                                   for (size_t i = start; i < end; ++i)
                                       it = i;
                               },
                               .number_of_user_bins = 500,
                               .maximum_fpr = 0.001,
                               .threads = 4,
                               .tmax = 64,
                               .disable_estimate_union = true,
                               .disable_rearrangement = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    auto agent = hibf.template counting_agent<uint32_t>();

    auto query = std::views::iota(1u, 251'000u); // 501^2 - 1
    auto & result = agent.bulk_count(query, 1u);

    for (size_t ub_id = 0u; ub_id < config.number_of_user_bins; ++ub_id)
    {
        // Hashes are consecutive integer values and are hence a bad input for the HIBF.
        // Only checking for minimum count.
        size_t const min_count = 2 * ub_id + 2; // size of [ squared(ub_id + 1u), ..., squared(ub_id + 2u) - 1u )
        ASSERT_GE(result[ub_id], min_count) << "ub_id: " << ub_id;
    }
}

TEST(hibf_test, serialisation)
{
    // range of range of sequences
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    seqan::hibf::config config{.input_fn =
                                   [&](size_t const num, seqan::hibf::insert_iterator it)
                               {
                                   for (auto const hash : hashes[num])
                                       it = hash;
                               },
                               .number_of_user_bins = 2};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

    seqan::hibf::test::test_serialisation(std::move(hibf));
}
