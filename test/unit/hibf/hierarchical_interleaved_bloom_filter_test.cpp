// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, TestInfo, TestPartResult, EXPECT_EQ

#include <cstddef>    // for size_t
#include <functional> // for function
#include <ranges>     // for _Iota, iota, views
#include <sstream>    // for basic_stringstream, stringstream
#include <vector>     // for vector, allocator

#include <hibf/config.hpp>                                // for insert_iterator, config
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/layout/layout.hpp>                         // for layout
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

TEST(hibf_test, set_number_of_user_bins)
{
    seqan::hibf::config config{.input_fn =
                                   [&](size_t const, seqan::hibf::insert_iterator it)
                               {
                                   it = 5u;
                               },
                               .number_of_user_bins = 73};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);

    hibf.set_number_of_user_bins();
    EXPECT_EQ(config.number_of_user_bins, hibf.number_of_user_bins);

    hibf.ibf_bin_to_user_bin_id.emplace_back(1u, 73);
    hibf.set_number_of_user_bins();
    EXPECT_EQ(config.number_of_user_bins + 1u, hibf.number_of_user_bins);
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

TEST(hibf_test, unevenly_sized_and_unique_user_bins)
{
    // std::pow(number, 2) does conversion do double and back to size_t.
    // If used more often, consider porting pow with integer overloads from
    // https://github.com/seqan/seqan3/blob/8b7d02cb8695369b7baeb4b3042ae7c864b67b8c/include/seqan3/utility/math.hpp
    auto squared = [](size_t const number)
    {
        return number * number;
    };

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

// #ifdef HIBF_HAS_SEQAN3

// #include <seqan3/alphabet/nucleotide/dna4.hpp>
// #include <seqan3/core/debug_stream.hpp>
// #include <seqan3/core/detail/all_view.hpp>
// #include <seqan3/io/sequence_file/all.hpp>
// #include <seqan3/test/cereal.hpp>
// #include <seqan3/utility/views/deep.hpp>

// TEST(hibf_seqan3_test, input_sequences_of_sequences)
// {
//     using namespace seqan::hibf::literals;

//     auto kmer_transformation = seqan::hibf::views::kmer_hash(seqan::hibf::ungapped{2u});

//     // range of range of sequences
//     std::vector<std::vector<std::vector<seqan::hibf::dna4>>> seqs{{"AAAGGGGGGC"_dna4}, {"TTTTTT"_dna4}};

//     seqan::hibf::config config
//     {
//         .input = seqs | seqan::hibf::views::deep{kmer_transformation},
//         .rearrange_user_bins = false
//     };

//     seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

//     auto agent = hibf.membership_agent();

//     std::vector<seqan::hibf::dna4> query{"AAGG"_dna4};
//     auto query_kmers = query | kmer_transformation;

//     auto result = agent.membership_for(query_kmers, 1);

//     seqan::hibf::debug_stream << result << std::endl;
// }

// TEST(hibf_seqan3_test, input_files)
// {
//     seqan::hibf::test::tmp_directory tmp{};
//     std::filesystem::path f1{tmp.path() / "f1.fa"};
//     std::filesystem::path f2{tmp.path() / "f2.fa"};

//     {
//         std::ofstream out{f1};
//         out << ">seq1\nAAAGGGGGGC\n";

//         std::ofstream out2{f2};
//         out2 << ">seq1\nTTTTTT\n";
//     }

//     auto transform = seqan::hibf::views::kmer_hash(seqan::hibf::ungapped{2u});

//     // range of range of sequences
//     std::vector<std::string> filenames{f1.string(), f2.string()};
//     auto file_range = filenames | std::views::transform([&transform](auto const & f)
//     {
//         auto record_transform = std::views::transform([&transform](auto && rec){ return rec.sequence() | transform; });
//         return seqan::hibf::detail::all(seqan::hibf::sequence_file_input{f}) | record_transform;
//     });

//     seqan::hibf::config config
//     {
//         .input = file_range,
//         .rearrange_user_bins = false
//     };

//     seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

//     auto agent = hibf.membership_agent();

//     using namespace seqan::hibf::literals;

//     std::vector<seqan::hibf::dna4> query{"AAGG"_dna4};

//     auto result = agent.membership_for(query | transform, 1);

//     seqan::hibf::debug_stream << result << std::endl; // prints [0] since query is found in user bin 0
// }

// #endif // HIBF_HAS_SEQAN3
