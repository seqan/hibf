// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <hibf/layout/compute_layout.hpp>
#include <hibf/misc/timer.hpp>
#include <hibf/sketch/compute_sketches.hpp>

TEST(compute_layout, dispatch)
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
    config.validate_and_set_defaults();

    std::vector<seqan::hibf::sketch::hyperloglog> sketches;
    std::vector<size_t> kmer_counts;

    seqan::hibf::sketch::compute_sketches(config, kmer_counts, sketches);

    auto layout1 = seqan::hibf::layout::compute_layout(config, kmer_counts, sketches);

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};

    auto layout2 =
        seqan::hibf::layout::compute_layout(config, kmer_counts, sketches, union_estimation_timer, rearrangement_timer);

    EXPECT_TRUE(layout1 == layout2);
}
