// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <hibf/test/expect_range_eq.hpp>
#include <hibf/test/seqan3.hpp>

TEST(hibf_test, test_specific_hash_values)
{
    // range of range of sequences
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    hibf::config config{.input_fn =
                            [&](size_t const num, hibf::insert_iterator it)
                        {
                            for (auto const hash : hashes[num])
                                it = hash;
                        },
                        .number_of_user_bins = 2,
                        .disable_rearrangement = true};

    hibf::hierarchical_interleaved_bloom_filter hibf{config};

    {
        std::vector<size_t> query{1, 2, 3, 4, 5};

        auto agent = hibf.membership_agent();
        auto result = agent.bulk_contains(query, 2);

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
//     using namespace hibf::literals;

//     auto kmer_transformation = hibf::views::kmer_hash(hibf::ungapped{2u});

//     // range of range of sequences
//     std::vector<std::vector<std::vector<hibf::dna4>>> seqs{{"AAAGGGGGGC"_dna4}, {"TTTTTT"_dna4}};

//     hibf::config config
//     {
//         .input = seqs | hibf::views::deep{kmer_transformation},
//         .rearrange_user_bins = false
//     };

//     hibf::hierarchical_interleaved_bloom_filter hibf{config};

//     auto agent = hibf.membership_agent();

//     std::vector<hibf::dna4> query{"AAGG"_dna4};
//     auto query_kmers = query | kmer_transformation;

//     auto result = agent.bulk_contains(query_kmers, 1);

//     hibf::debug_stream << result << std::endl;
// }

// TEST(hibf_seqan3_test, input_files)
// {
//     hibf::test::tmp_directory tmp{};
//     std::filesystem::path f1{tmp.path() / "f1.fa"};
//     std::filesystem::path f2{tmp.path() / "f2.fa"};

//     {
//         std::ofstream out{f1};
//         out << ">seq1\nAAAGGGGGGC\n";

//         std::ofstream out2{f2};
//         out2 << ">seq1\nTTTTTT\n";
//     }

//     auto transform = hibf::views::kmer_hash(hibf::ungapped{2u});

//     // range of range of sequences
//     std::vector<std::string> filenames{f1.string(), f2.string()};
//     auto file_range = filenames | std::views::transform([&transform](auto const & f)
//     {
//         auto record_transform = std::views::transform([&transform](auto && rec){ return rec.sequence() | transform; });
//         return hibf::detail::all(hibf::sequence_file_input{f}) | record_transform;
//     });

//     hibf::config config
//     {
//         .input = file_range,
//         .rearrange_user_bins = false
//     };

//     hibf::hierarchical_interleaved_bloom_filter hibf{config};

//     auto agent = hibf.membership_agent();

//     using namespace hibf::literals;

//     std::vector<hibf::dna4> query{"AAGG"_dna4};

//     auto result = agent.bulk_contains(query | transform, 1);

//     hibf::debug_stream << result << std::endl; // prints [0] since query is found in user bin 0
// }

// #endif // HIBF_HAS_SEQAN3
