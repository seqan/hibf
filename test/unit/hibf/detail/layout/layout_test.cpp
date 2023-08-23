#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <hibf/detail/layout/layout.hpp> // for layout, operator<<

TEST(layout_test, printing_max_bins)
{
    std::stringstream ss{};

    hibf::layout::layout layout;

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

    hibf::layout::layout layout;

    layout.user_bins.emplace_back(7, std::vector<size_t>{}, 1, 0);
    layout.user_bins.emplace_back(4, std::vector<size_t>{1}, 22, 0);
    layout.user_bins.emplace_back(5, std::vector<size_t>{1, 2, 3, 4}, 21, 22);

    for (auto const & ub : layout.user_bins)
        ss << ub << "\n";

    std::string expected = R"ub(7	0	1
4	1;0	1;22
5	1;2;3;4;22	1;1;1;1;21
)ub";

    EXPECT_EQ(ss.str(), expected);
}
