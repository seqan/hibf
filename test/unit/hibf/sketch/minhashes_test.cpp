// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for AssertionResult, Message, Test, TestPartResult, EXPECT_FALSE, TestInfo

#include <algorithm> // for __fn, equal, make_heap
#include <cinttypes> // for uint64_t
#include <cstddef>   // for size_t
#include <ranges>    // for __fn, iota, views
#include <span>      // for span
#include <utility>   // for move
#include <vector>    // for vector

#include <hibf/misc/iota_vector.hpp>     // for iota_vector
#include <hibf/sketch/minhashes.hpp>     // for minhashes
#include <hibf/test/expect_range_eq.hpp> // for expect_range_eq, EXPECT_RANGE_EQ

TEST(minhashes_test, ctor)
{
    seqan::hibf::sketch::minhashes default_ctor{};
    EXPECT_TRUE(default_ctor.table.empty());
    EXPECT_FALSE(default_ctor.is_valid());

    seqan::hibf::sketch::minhashes copy_ctor{default_ctor};
    EXPECT_TRUE(copy_ctor.table.empty());
    EXPECT_FALSE(copy_ctor.is_valid());

    seqan::hibf::sketch::minhashes copy_assignment{};
    copy_assignment = default_ctor;
    EXPECT_TRUE(copy_assignment.table.empty());
    EXPECT_FALSE(copy_assignment.is_valid());

    seqan::hibf::sketch::minhashes move_ctor{std::move(default_ctor)};
    EXPECT_TRUE(move_ctor.table.empty());
    EXPECT_FALSE(move_ctor.is_valid());

    seqan::hibf::sketch::minhashes move_assignment{};
    move_assignment = std::move(copy_ctor);
    EXPECT_TRUE(move_assignment.table.empty());
    EXPECT_FALSE(move_assignment.is_valid());
}

TEST(minhashes_test, push_to_heap_if_smaller)
{
    std::vector<uint64_t> heap{3, 4, 5, 6, 7, 8};
    std::ranges::make_heap(heap);

    // No change because 10 is larger than all values in the heap
    std::vector<uint64_t> const original_heap{heap};
    seqan::hibf::sketch::minhashes::push_to_heap_if_smaller(10, heap);
    EXPECT_RANGE_EQ(original_heap, heap);

    // Heap changes
    seqan::hibf::sketch::minhashes::push_to_heap_if_smaller(0, heap);
    EXPECT_FALSE(std::ranges::equal(original_heap, heap));
    EXPECT_EQ(heap[0], 7u); // largest element is now 7 (not 8 anymore)
}

TEST(minhashes_test, ctor_sorted_list)
{
    std::vector<uint64_t> const heap = seqan::hibf::iota_vector<uint64_t>(1000);

    seqan::hibf::sketch::minhashes sketches{heap};

    // check for correct construction
    EXPECT_FALSE(sketches.table.empty());
    EXPECT_TRUE(sketches.is_valid());

    ASSERT_EQ(sketches.table.size(), seqan::hibf::sketch::minhashes::num_sketches);

    size_t sid{};
    for (auto & sketch : sketches.table)
    {
        ASSERT_EQ(sketch.size(), seqan::hibf::sketch::minhashes::sketch_size) << "sid: " << sid;
        EXPECT_RANGE_EQ(sketch, std::views::iota(0u, seqan::hibf::sketch::minhashes::sketch_size)) << "sid: " << sid;
        ++sid;
    }
}

TEST(minhashes_test, fill_incomplete_sketches)
{
    std::vector<uint64_t> const heap = seqan::hibf::iota_vector<uint64_t>(10);

    seqan::hibf::sketch::minhashes sketches{heap}; // there will be too few values
    EXPECT_FALSE(sketches.table.empty());
    ASSERT_FALSE(sketches.is_valid());

    // recreate bigger heap
    std::vector<uint64_t> bigger_heap = seqan::hibf::iota_vector<uint64_t>(1000);
    std::span<uint64_t> const bigger_heap_non_redundant(bigger_heap.begin() + heap.size(), bigger_heap.end());
    sketches.fill_incomplete_sketches(bigger_heap_non_redundant);

    // check if filling up worked correctly
    ASSERT_TRUE(sketches.is_valid());

    ASSERT_EQ(sketches.table.size(), seqan::hibf::sketch::minhashes::num_sketches);

    size_t sid{};
    for (auto & sketch : sketches.table)
    {
        ASSERT_EQ(sketch.size(), seqan::hibf::sketch::minhashes::sketch_size) << "sid: " << sid;
        EXPECT_RANGE_EQ(sketch, std::views::iota(0u, seqan::hibf::sketch::minhashes::sketch_size)) << "sid: " << sid;
        ++sid;
    }
}
