// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestPartResult, Test, TestInfo, EXPECT_NO_THROW, TEST

#include <cstddef>   // for size_t
#include <stdexcept> // for invalid_argument
#include <string>    // for basic_string
#include <vector>    // for vector

#include <hibf/layout/data_store.hpp>     // for data_store
#include <hibf/layout/layout.hpp>         // for layout
#include <hibf/sketch/hyperloglog.hpp>    // for hyperloglog
#include <hibf/test/expect_throw_msg.hpp> // for EXPECT_THROW_MSG

TEST(data_store_test, validate)
{
    seqan::hibf::layout::layout layout{};
    std::vector<size_t> kmer_counts(3);
    std::vector<seqan::hibf::sketch::hyperloglog> sketches(3);

    // hibf_layout must not be nullptr
    {
        seqan::hibf::layout::data_store store{};
        EXPECT_THROW_MSG(store.validate(),
                         std::invalid_argument,
                         "[HIBF ERROR] data_store::hibf_layout must not be nullptr.");
    }

    // kmer_counts must not be nullptr
    {
        seqan::hibf::layout::data_store store{.hibf_layout = &layout};
        EXPECT_THROW_MSG(store.validate(),
                         std::invalid_argument,
                         "[HIBF ERROR] data_store::kmer_counts must not be nullptr.");
    }

    // kmer_counts and sketches must have the same size
    {
        std::vector<seqan::hibf::sketch::hyperloglog> wrong_sketches(2);

        seqan::hibf::layout::data_store store{.hibf_layout = &layout,
                                              .kmer_counts = &kmer_counts,
                                              .sketches = &wrong_sketches};
        EXPECT_THROW_MSG(store.validate(),
                         std::invalid_argument,
                         "[HIBF ERROR] data_store::kmer_counts and data_store::sketches must have the same size.");
    }

    // kmer_counts size must be greater than positions size
    {
        seqan::hibf::layout::data_store store{.hibf_layout = &layout,
                                              .kmer_counts = &kmer_counts,
                                              .sketches = &sketches,
                                              .positions = {1, 2, 3, 4}};

        EXPECT_THROW_MSG(
            store.validate(),
            std::invalid_argument,
            "[HIBF ERROR] data_store::kmer_counts.size() must not be smaller than data_store::positions.size().");
    }

    // fpr_correction must not be empty
    {
        seqan::hibf::layout::data_store store{.hibf_layout = &layout,
                                              .kmer_counts = &kmer_counts,
                                              .sketches = &sketches};

        EXPECT_THROW_MSG(store.validate(),
                         std::invalid_argument,
                         "[HIBF ERROR] data_store::fpr_correction must not be empty.");
    }

    // relaxed_fpr_correction must be in (0.0,1.0]
    {
        seqan::hibf::layout::data_store store{.hibf_layout = &layout,
                                              .kmer_counts = &kmer_counts,
                                              .sketches = &sketches,
                                              .fpr_correction = {1.0, 2.0, 3.0}};

        EXPECT_THROW_MSG(store.validate(),
                         std::invalid_argument,
                         "[HIBF ERROR] data_store::relaxed_fpr_correction must be in (0.0,1.0].");

        store.relaxed_fpr_correction = 1.01;
        EXPECT_THROW_MSG(store.validate(),
                         std::invalid_argument,
                         "[HIBF ERROR] data_store::relaxed_fpr_correction must be in (0.0,1.0].");

        store.relaxed_fpr_correction = 0.5;
        EXPECT_NO_THROW(store.validate());
    }
}
