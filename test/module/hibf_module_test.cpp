// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Checks that the public API is usable through `import seqan.hibf;`.
 * \details
 * Every textual `#include` must appear before the `import`: GCC 16 rejects the reverse order (it sees the standard
 * library declarations of the global module fragment a second time). Clang accepts both orders.
 */

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <sstream>
#include <string>
#include <vector>

#include <cereal/archives/binary.hpp>

import seqan.hibf;

TEST(module_test, version)
{
    EXPECT_EQ(seqan::hibf::hibf_version,
              seqan::hibf::hibf_version_major * 10000 + seqan::hibf::hibf_version_minor * 100
                  + seqan::hibf::hibf_version_patch);
    EXPECT_NE(std::string{seqan::hibf::hibf_version_cstring}, std::string{});
}

TEST(module_test, free_functions)
{
    EXPECT_EQ(seqan::hibf::next_multiple_of_64(1u), 64u);
    EXPECT_EQ(seqan::hibf::divide_and_ceil(5u, 2u), 3u);
    EXPECT_EQ(seqan::hibf::add_empty_bins(64u, 0.5), 128u);
    EXPECT_EQ(seqan::hibf::subtract_empty_bins(128u, 0.5), 64u);

    // constexpr evaluation of an imported function template
    static_assert(seqan::hibf::iota_vector<uint16_t>(3u).size() == 3u);
    EXPECT_EQ(seqan::hibf::iota_vector(3u), (std::vector<size_t>{0u, 1u, 2u}));
}

TEST(module_test, constants)
{
    EXPECT_EQ(seqan::hibf::prefix::meta_header, "@");
    EXPECT_TRUE(seqan::hibf::prefix::meta_hibf_config_start.starts_with(seqan::hibf::prefix::meta_header));
    EXPECT_GT(seqan::hibf::bin_kind::merged, seqan::hibf::bin_kind::deleted);
}

TEST(module_test, interleaved_bloom_filter)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{64u},
                                              seqan::hibf::bin_size{1024u},
                                              seqan::hibf::hash_function_count{2u}};

    EXPECT_EQ(ibf.bin_count(), 64u);
    EXPECT_EQ(ibf.bin_size(), 1024u);
    EXPECT_EQ(ibf.hash_function_count(), 2u);

    ibf.emplace(126u, seqan::hibf::bin_index{0u});
    ibf.emplace(712u, seqan::hibf::bin_index{0u});
    ibf.emplace(126u, seqan::hibf::bin_index{3u});

    auto containment_agent = ibf.containment_agent();
    seqan::hibf::bit_vector const & containment = containment_agent.bulk_contains(126u);
    EXPECT_TRUE(containment[0u]);
    EXPECT_TRUE(containment[3u]);
    EXPECT_FALSE(containment[1u]);

    auto counting_agent = ibf.counting_agent();
    seqan::hibf::counting_vector<uint16_t> const & counts = counting_agent.bulk_count(std::vector<size_t>{126u, 712u});
    EXPECT_EQ(counts[0u], 2u);

    auto membership_agent = ibf.membership_agent();
    EXPECT_EQ(membership_agent.membership_for(std::vector<size_t>{126u, 712u}, 2u), (std::vector<uint64_t>{0u}));
}

TEST(module_test, hierarchical_interleaved_bloom_filter)
{
    std::vector<std::vector<size_t>> const hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    // seqan::hibf::insert_iterator is part of the public API and is used through the imported config.
    auto input = [&hashes](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        for (size_t const hash : hashes[user_bin_id])
            it = hash;
    };

    seqan::hibf::config config{.input_fn = input, .number_of_user_bins = 2u};
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

    auto agent = hibf.membership_agent();
    std::vector<uint64_t> const & result = agent.membership_for(std::vector<size_t>{8u, 9u, 10u}, 3u);
    EXPECT_EQ(result, (std::vector<uint64_t>{0u}));

    // seqan::hibf::print is a function object, i.e. an exported variable, not a function.
    std::ostringstream stream{};
    seqan::hibf::print(result, stream);
    EXPECT_EQ(stream.str(), "[0]\n");
}

TEST(module_test, serialisation)
{
    seqan::hibf::interleaved_bloom_filter original{seqan::hibf::bin_count{64u}, seqan::hibf::bin_size{1024u}};
    original.emplace(126u, seqan::hibf::bin_index{0u});

    std::stringstream stream{};
    {
        cereal::BinaryOutputArchive archive{stream};
        archive(original);
    }

    seqan::hibf::interleaved_bloom_filter restored{};
    {
        cereal::BinaryInputArchive archive{stream};
        archive(restored);
    }

    EXPECT_TRUE(original == restored);
}

TEST(module_test, path_serialisation)
{
    // hibf/cereal/path.hpp adds `save`/`load` overloads in namespace cereal. cereal finds them by ADL, so they have to
    // be visible to importers, not merely reachable. See the `export namespace cereal` block in src/seqan.hibf.cppm.
    std::filesystem::path const original{"/some/random/path.txt"};

    std::stringstream stream{};
    {
        cereal::BinaryOutputArchive archive{stream};
        archive(original);
    }

    std::filesystem::path restored{};
    {
        cereal::BinaryInputArchive archive{stream};
        archive(restored);
    }

    EXPECT_EQ(original, restored);
}

TEST(module_test, layout)
{
    // Hidden friends of imported types must still be found by ADL.
    seqan::hibf::layout::layout::user_bin const user_bin{.previous_TB_indices = {},
                                                         .storage_TB_id = 0u,
                                                         .number_of_technical_bins = 1u,
                                                         .idx = 0u};

    std::ostringstream stream{};
    stream << user_bin;
    EXPECT_EQ(stream.str(), "0\t0\t1");

    seqan::hibf::layout::layout layout{};
    layout.user_bins.push_back(user_bin);
    EXPECT_EQ(layout.user_bins.size(), 1u);

    seqan::hibf::config const config{.number_of_user_bins = 1u};
    EXPECT_TRUE(layout.validate(config));
    layout.top_level_max_bin_id = 1u;
    std::ostringstream diagnostic_stream{};
    std::string formatted{};
    EXPECT_FALSE(layout.validate(config,
                                 [&](seqan::hibf::layout::layout::diagnostic const & diagnostic)
                                 {
                                     diagnostic_stream << diagnostic;
                                     formatted = std::format("{}", diagnostic);
                                 }));
    EXPECT_TRUE(diagnostic_stream.str().starts_with("[HIBF LAYOUT ERROR] "));
    EXPECT_EQ(formatted, diagnostic_stream.str());

    std::vector<double> const correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2u, .t_max = 4u});
    EXPECT_EQ(correction.size(), seqan::hibf::next_multiple_of_64(4u) + 1u);
    EXPECT_DOUBLE_EQ(correction[1u], 1.0);
    EXPECT_GT(correction[2u], 1.0);
}

TEST(module_test, sketch)
{
    seqan::hibf::sketch::hyperloglog sketch{12u};
    EXPECT_EQ(sketch.data_size(), 1u << 12u);
    EXPECT_DOUBLE_EQ(sketch.estimate(), 0.0);

    for (size_t i = 0u; i < 1000u; ++i)
        sketch.add(i);

    // HyperLogLog only approximates; this test is about reaching the API through the module, not about its accuracy.
    double const estimate = sketch.estimate();
    EXPECT_GT(estimate, 500.0);
    EXPECT_LT(estimate, 2000.0);

    seqan::hibf::sketch::hyperloglog other{12u};
    for (size_t i = 500u; i < 1500u; ++i)
        other.add(i);

    // `sketch` holds 1000 and `other` holds 1000 distinct values, overlapping in 500. The union is larger than either.
    EXPECT_GT(sketch.merge_and_estimate(other), estimate);
}

TEST(module_test, timer)
{
    seqan::hibf::serial_timer timer{};
    timer.start();
    timer.stop();
    EXPECT_GE(timer.in_seconds(), 0.0);

    seqan::hibf::concurrent_timer concurrent{};
    concurrent += timer;
    EXPECT_GE(concurrent.in_seconds(), 0.0);
    EXPECT_DOUBLE_EQ(concurrent.avg_in_seconds(), concurrent.in_seconds());
}

TEST(module_test, concepts)
{
    static_assert(seqan::hibf::cereal_output_archive<cereal::BinaryOutputArchive>);
    static_assert(seqan::hibf::cereal_input_archive<cereal::BinaryInputArchive>);
    static_assert(seqan::hibf::cereal_archive<cereal::BinaryOutputArchive>);
    static_assert(!seqan::hibf::cereal_text_archive<cereal::BinaryOutputArchive>);
}
