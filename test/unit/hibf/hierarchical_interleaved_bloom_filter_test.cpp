// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, Message, TEST, TestPartResult

#include <cstddef>    // for size_t
#include <functional> // for function
#include <sstream>    // for basic_stringstream, stringstream
#include <vector>     // for vector, allocator

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/layout/layout.hpp>                         // for layout
#include <hibf/test/expect_range_eq.hpp>                  // for expect_range_eq, EXPECT_RANGE_EQ

TEST(hibf_test, test_specific_hash_values)
{
    // range of range of sequences
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    seqan::hibf::config config{.input_fn =
                                   [&](size_t const num, seqan::hibf::insert_iterator it)
                               {
                                   for (auto const hash : hashes[num])
                                       it = hash;
                               },
                               .number_of_user_bins = 2,
                               .disable_rearrangement = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

    {
        std::vector<size_t> query{1, 2, 3, 4, 5};

        auto agent = hibf.membership_agent();
        auto result = agent.membership_for(query, 2);

        EXPECT_RANGE_EQ(result, (std::vector<size_t>{0u, 1u}));
    }
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
                             "@        \"maximum_false_positive_rate\": 0.05,\n"
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

    {
        std::vector<size_t> query{1, 2, 3, 4, 5};

        auto agent = hibf.membership_agent();
        auto result = agent.membership_for(query, 2);

        EXPECT_RANGE_EQ(result, (std::vector<size_t>{0u, 1u}));
    }
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
