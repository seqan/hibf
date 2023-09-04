#include <gtest/gtest.h> // for Test, EXPECT_EQ, Message, TestPartResult, TestInfo, TEST

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, basic_stringstream, basic_ostream, stringstream
#include <string>      // for allocator, char_traits, string
#include <string_view> // for operator<<
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

    std::string expected = R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)layout_file";

    EXPECT_EQ(ss.str(), expected);
}

TEST(layout_test, read_from)
{
    std::stringstream ss{R"layout_file(#TOP_LEVEL_IBF fullest_technical_bin_idx:111
#LOWER_LEVEL_IBF_0 fullest_technical_bin_idx:0
#LOWER_LEVEL_IBF_2 fullest_technical_bin_idx:2
#LOWER_LEVEL_IBF_1;2;3;4 fullest_technical_bin_idx:22
#USER_BIN_IDX	TECHNICAL_BIN_INDICES	NUMBER_OF_TECHNICAL_BINS
7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)layout_file"};

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
