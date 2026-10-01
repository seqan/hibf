// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, AssertionResult, TestPartResult, CmpHelperEQ, CmpHelperEQFa...

#include <cstddef>     // for size_t
#include <sstream>     // for basic_stringstream, operator<<, stringstream, basic_ios, basic_iostream
#include <string>      // for char_traits, allocator, basic_string, string
#include <string_view> // for basic_string_view
#include <vector>      // for vector

#include <hibf/layout/layout.hpp> // for layout, operator<<

TEST(layout_test, printing_max_bins)
{
    std::stringstream ss{};

    seqan::hibf::layout::layout layout;

    layout.max_bins.emplace_back(std::vector<size_t>{}, 0);
    layout.max_bins.emplace_back(std::vector<size_t>{2}, 2);
    layout.max_bins.emplace_back(std::vector<size_t>{1, 2, 3, 4}, 22);

    for (auto const & mb : layout.max_bins)
        ss << mb << "\n";

    std::string expected = R"mb(#LOWER_LEVEL_IBF_ fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
)mb";

    EXPECT_EQ(ss.str(), expected);
}

TEST(layout_test, printing_user_bins)
{
    std::stringstream ss{};

    seqan::hibf::layout::layout layout;

    layout.user_bins.emplace_back(std::vector<size_t>{}, 0, 1, 7);
    layout.user_bins.emplace_back(std::vector<size_t>{1}, 0, 22, 4);
    layout.user_bins.emplace_back(std::vector<size_t>{1, 2, 3, 4}, 22, 21, 5);

    for (auto const & ub : layout.user_bins)
        ss << ub << "\n";

    std::string expected = R"ub(7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)ub";

    EXPECT_EQ(ss.str(), expected);
}

static std::string const layout_file{
    R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)layout_file"};

TEST(layout_test, write_to)
{
    std::stringstream ss{};

    seqan::hibf::layout::layout layout;

    layout.top_level_max_bin_id = 111;
    layout.max_bins.emplace_back(std::vector<size_t>{0}, 0);
    layout.max_bins.emplace_back(std::vector<size_t>{2}, 2);
    layout.max_bins.emplace_back(std::vector<size_t>{1, 2, 3, 4}, 22);
    layout.user_bins.emplace_back(std::vector<size_t>{}, 0, 1, 7);
    layout.user_bins.emplace_back(std::vector<size_t>{1}, 0, 22, 4);
    layout.user_bins.emplace_back(std::vector<size_t>{1, 2, 3, 4}, 22, 21, 5);

    layout.write_to(ss);

    EXPECT_EQ(ss.str(), layout_file);
}

TEST(layout_test, read_from)
{
    std::stringstream ss{layout_file};

    seqan::hibf::layout::layout layout;
    layout.read_from(ss);

    EXPECT_EQ(layout.top_level_max_bin_id, 111);
    EXPECT_EQ(layout.max_bins[0], (seqan::hibf::layout::layout::max_bin{{0}, 0}));
    EXPECT_EQ(layout.max_bins[1], (seqan::hibf::layout::layout::max_bin{{2}, 2}));
    EXPECT_EQ(layout.max_bins[2], (seqan::hibf::layout::layout::max_bin{{1, 2, 3, 4}, 22}));
    EXPECT_EQ(layout.user_bins[0], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 0, 1, 7}));
    EXPECT_EQ(layout.user_bins[1], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 0, 22, 4}));
    EXPECT_EQ(layout.user_bins[2], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1, 2, 3, 4}, 22, 21, 5}));
}

TEST(layout_test, clear)
{
    std::stringstream ss{layout_file};

    seqan::hibf::layout::layout layout;
    layout.read_from(ss);

    ASSERT_NE(layout.top_level_max_bin_id, 0);
    ASSERT_FALSE(layout.max_bins.empty());
    ASSERT_FALSE(layout.user_bins.empty());

    layout.clear();

    EXPECT_EQ(layout.top_level_max_bin_id, 0);
    EXPECT_TRUE(layout.max_bins.empty());
    EXPECT_TRUE(layout.user_bins.empty());
}

TEST(layout_test, read_from_partitioned_layout)
{
    // layout consists of three partitions, written one after the other
    std::stringstream ss{R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:1
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)layout_file"};

    for (size_t i = 0; i < 3; ++i)
    {
        seqan::hibf::layout::layout layout;
        layout.read_from(ss);

        EXPECT_EQ(layout.top_level_max_bin_id, 111);
        EXPECT_EQ(layout.max_bins[0], (seqan::hibf::layout::layout::max_bin{{0}, i}));
        EXPECT_EQ(layout.max_bins[1], (seqan::hibf::layout::layout::max_bin{{2}, 2}));
        EXPECT_EQ(layout.max_bins[2], (seqan::hibf::layout::layout::max_bin{{1, 2, 3, 4}, 22}));
        EXPECT_EQ(layout.user_bins[0], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 0, 1, 7}));
        EXPECT_EQ(layout.user_bins[1], (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 0, 22, 4}));
        EXPECT_EQ(layout.user_bins[2],
                  (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1, 2, 3, 4}, 22, 21, 5}));
    }
}

namespace
{

using layout_t = seqan::hibf::layout::layout;

// Root-IBF: UB 0 is split into TBs 0-1, TB 2 is merged, UB 5 is split into TBs 3-63.
// IBF 2: UB 1 is split into TBs 0-1, TB 2 is merged, UB 2 is split into TBs 3-63.
// IBF 2;2: UB 3 in TB 0, UB 4 is split into TBs 1-63.
// Each IBF uses all of its 64 technical bins. The user bins are deliberately not sorted by index.
layout_t valid_layout()
{
    layout_t layout{};
    layout.top_level_max_bin_id = 2;
    layout.max_bins.emplace_back(std::vector<size_t>{2}, 0);
    layout.max_bins.emplace_back(std::vector<size_t>{2, 2}, 1);
    layout.user_bins.emplace_back(std::vector<size_t>{}, 3, 61, 5);
    layout.user_bins.emplace_back(std::vector<size_t>{}, 0, 2, 0);
    layout.user_bins.emplace_back(std::vector<size_t>{2}, 0, 2, 1);
    layout.user_bins.emplace_back(std::vector<size_t>{2}, 3, 61, 2);
    layout.user_bins.emplace_back(std::vector<size_t>{2, 2}, 0, 1, 3);
    layout.user_bins.emplace_back(std::vector<size_t>{2, 2}, 1, 63, 4);
    return layout;
}

} // namespace

TEST(layout_test, number_of_levels)
{
    layout_t layout = valid_layout();
    EXPECT_EQ(layout.number_of_levels(), 3u);

    layout.user_bins.resize(2u); // only the top-level user bins
    EXPECT_EQ(layout.number_of_levels(), 1u);

    layout.clear();
    EXPECT_EQ(layout.number_of_levels(), 0u);
}
