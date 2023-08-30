#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <hibf/layout/graph.hpp>         // for layout, operator<<
#include <hibf/test/expect_range_eq.hpp> // for expect_range_eq, EXPECT_RANGE_EQ

TEST(layout_test, printing_max_bins)
{
    // prepare layout
    seqan::hibf::layout::layout hibf_layout;
    // set top level max bin id
    hibf_layout.top_level_max_bin_id = 0;
    // set max bin id per lower level ibf
    hibf_layout.max_bins.emplace_back(std::vector<size_t>{0}, 1);
    hibf_layout.max_bins.emplace_back(std::vector<size_t>{1}, 26);
    hibf_layout.max_bins.emplace_back(std::vector<size_t>{0, 0}, 4);
    hibf_layout.max_bins.emplace_back(std::vector<size_t>{0, 1}, 34);
    hibf_layout.max_bins.emplace_back(std::vector<size_t>{0, 0, 0}, 30);
    // isnert (previous_TB_indices, storage_TB_id, number_of_technical_bins, idx)
    hibf_layout.user_bins.emplace_back(15, std::vector<size_t>{0, 0, 0}, 30, 0);
    hibf_layout.user_bins.emplace_back(16, std::vector<size_t>{0, 0, 0}, 11, 30);
    hibf_layout.user_bins.emplace_back(17, std::vector<size_t>{0, 0, 0}, 11, 41);
    hibf_layout.user_bins.emplace_back(18, std::vector<size_t>{0, 0, 0}, 6, 52);
    hibf_layout.user_bins.emplace_back(19, std::vector<size_t>{0, 0, 0}, 6, 58);
    hibf_layout.user_bins.emplace_back(14, std::vector<size_t>{0, 0}, 1, 1);
    hibf_layout.user_bins.emplace_back(13, std::vector<size_t>{0, 0}, 1, 2);
    hibf_layout.user_bins.emplace_back(12, std::vector<size_t>{0, 0}, 1, 3);
    hibf_layout.user_bins.emplace_back(11, std::vector<size_t>{0, 0}, 1, 4);
    hibf_layout.user_bins.emplace_back(8, std::vector<size_t>{0, 1}, 34, 0);
    hibf_layout.user_bins.emplace_back(9, std::vector<size_t>{0, 1}, 15, 34);
    hibf_layout.user_bins.emplace_back(10, std::vector<size_t>{0, 1}, 15, 49);
    hibf_layout.user_bins.emplace_back(7, std::vector<size_t>{0}, 1, 2);
    hibf_layout.user_bins.emplace_back(6, std::vector<size_t>{0}, 1, 3);
    hibf_layout.user_bins.emplace_back(5, std::vector<size_t>{0}, 1, 4);
    hibf_layout.user_bins.emplace_back(2, std::vector<size_t>{1}, 26, 0);
    hibf_layout.user_bins.emplace_back(3, std::vector<size_t>{1}, 19, 26);
    hibf_layout.user_bins.emplace_back(4, std::vector<size_t>{1}, 19, 45);
    hibf_layout.user_bins.emplace_back(1, std::vector<size_t>{}, 1, 2);
    hibf_layout.user_bins.emplace_back(0, std::vector<size_t>{}, 2, 3);

    seqan::hibf::layout::graph hibf_graph{hibf_layout};

    // root node
    EXPECT_EQ(hibf_graph.root.parent_bin_index, 0);
    EXPECT_EQ(hibf_graph.root.max_bin_index, 0);
    EXPECT_EQ(hibf_graph.root.number_of_technical_bins, 5);
    EXPECT_EQ(hibf_graph.root.favourite_child_idx, 0);
    ASSERT_EQ(hibf_graph.root.remaining_records.size(), 2);
    EXPECT_EQ(hibf_graph.root.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{1, std::vector<size_t>{}, 1, 2}));
    EXPECT_EQ(hibf_graph.root.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{0, std::vector<size_t>{}, 2, 3}));

    // tree structure:
    ASSERT_EQ(hibf_graph.root.children.size(), 2);
    auto const & root_child0 = hibf_graph.root.children[0];
    auto const & root_child1 = hibf_graph.root.children[1];
    ASSERT_EQ(root_child0.children.size(), 2);
    auto const & c0_child0 = root_child0.children[0];
    auto const & c0_child1 = root_child0.children[1];
    ASSERT_EQ(c0_child0.children.size(), 1);
    auto const & c0_c0_child0 = c0_child0.children[0];

    // now test children sperarately

    // child 0 of root
    EXPECT_EQ(root_child0.parent_bin_index, 0);
    EXPECT_EQ(root_child0.max_bin_index, 1);
    EXPECT_EQ(root_child0.number_of_technical_bins, 5);
    EXPECT_EQ(root_child0.favourite_child_idx, 1);
    ASSERT_EQ(root_child0.remaining_records.size(), 3);
    EXPECT_EQ(root_child0.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{7, std::vector<size_t>{0}, 1, 2}));
    EXPECT_EQ(root_child0.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{6, std::vector<size_t>{0}, 1, 3}));
    EXPECT_EQ(root_child0.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{5, std::vector<size_t>{0}, 1, 4}));

    // child 1 of root
    EXPECT_EQ(root_child1.parent_bin_index, 1);
    EXPECT_EQ(root_child1.max_bin_index, 26);
    EXPECT_EQ(root_child1.number_of_technical_bins, 64);
    EXPECT_EQ(root_child1.favourite_child_idx, -1);
    ASSERT_EQ(root_child1.remaining_records.size(), 3);
    EXPECT_EQ(root_child1.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{3, std::vector<size_t>{1}, 19, 26}));
    EXPECT_EQ(root_child1.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{2, std::vector<size_t>{1}, 26, 0}));
    EXPECT_EQ(root_child1.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{4, std::vector<size_t>{1}, 19, 45}));

    // child 0 of child 0
    EXPECT_EQ(c0_child0.parent_bin_index, 0);
    EXPECT_EQ(c0_child0.max_bin_index, 4);
    EXPECT_EQ(c0_child0.number_of_technical_bins, 5);
    EXPECT_EQ(c0_child0.favourite_child_idx, -1);
    ASSERT_EQ(c0_child0.remaining_records.size(), 4);
    EXPECT_EQ(c0_child0.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{11, std::vector<size_t>{0, 0}, 1, 4}));
    EXPECT_EQ(c0_child0.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{14, std::vector<size_t>{0, 0}, 1, 1}));
    EXPECT_EQ(c0_child0.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{13, std::vector<size_t>{0, 0}, 1, 2}));
    EXPECT_EQ(c0_child0.remaining_records[3],
              (seqan::hibf::layout::layout::user_bin{12, std::vector<size_t>{0, 0}, 1, 3}));

    // child 1 of child 0
    EXPECT_EQ(c0_child1.parent_bin_index, 1);
    EXPECT_EQ(c0_child1.max_bin_index, 34);
    EXPECT_EQ(c0_child1.number_of_technical_bins, 64);
    EXPECT_EQ(c0_child1.favourite_child_idx, -1);
    ASSERT_EQ(c0_child1.remaining_records.size(), 3);
    EXPECT_EQ(c0_child1.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{9, std::vector<size_t>{0, 1}, 15, 34}));
    EXPECT_EQ(c0_child1.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{8, std::vector<size_t>{0, 1}, 34, 0}));
    EXPECT_EQ(c0_child1.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{10, std::vector<size_t>{0, 1}, 15, 49}));

    // child 0 of child 0 of child 0
    EXPECT_EQ(c0_c0_child0.parent_bin_index, 0);
    EXPECT_EQ(c0_c0_child0.max_bin_index, 30);
    EXPECT_EQ(c0_c0_child0.number_of_technical_bins, 64);
    EXPECT_EQ(c0_c0_child0.favourite_child_idx, -1);
    ASSERT_EQ(c0_c0_child0.remaining_records.size(), 5);
    EXPECT_EQ(c0_c0_child0.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{16, std::vector<size_t>{0, 0, 0}, 11, 30}));
    EXPECT_EQ(c0_c0_child0.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{15, std::vector<size_t>{0, 0, 0}, 30, 0}));
    EXPECT_EQ(c0_c0_child0.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{17, std::vector<size_t>{0, 0, 0}, 11, 41}));
    EXPECT_EQ(c0_c0_child0.remaining_records[3],
              (seqan::hibf::layout::layout::user_bin{18, std::vector<size_t>{0, 0, 0}, 6, 52}));
    EXPECT_EQ(c0_c0_child0.remaining_records[4],
              (seqan::hibf::layout::layout::user_bin{19, std::vector<size_t>{0, 0, 0}, 6, 58}));
}
