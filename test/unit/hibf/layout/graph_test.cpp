#include <gtest/gtest.h> // for Message, TestPartResult, EXPECT_EQ, Test, ASSERT_EQ, TestInfo, TEST

#include <cstddef>  // for size_t
#include <optional> // for optional
#include <vector>   // for vector, allocator

#include <hibf/layout/graph.hpp>  // for graph
#include <hibf/layout/layout.hpp> // for layout

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
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0, 0}, 0, 30, 15);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0, 0}, 30, 11, 16);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0, 0}, 41, 11, 17);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0, 0}, 52, 6, 18);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0, 0}, 58, 6, 19);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0}, 1, 1, 14);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0}, 2, 1, 13);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0}, 3, 1, 12);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 0}, 4, 1, 11);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 1}, 0, 34, 8);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 1}, 34, 15, 9);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0, 1}, 49, 15, 10);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0}, 2, 1, 7);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0}, 3, 1, 6);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{0}, 4, 1, 5);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{1}, 0, 26, 2);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{1}, 26, 19, 3);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{1}, 45, 19, 4);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{}, 2, 1, 1);
    hibf_layout.user_bins.emplace_back(std::vector<size_t>{}, 3, 2, 0);

    seqan::hibf::layout::graph hibf_graph{hibf_layout};

    // root node
    EXPECT_EQ(hibf_graph.root.parent_bin_index, 0);
    EXPECT_EQ(hibf_graph.root.max_bin_index, 0);
    EXPECT_EQ(hibf_graph.root.number_of_technical_bins, 5);
    EXPECT_EQ(hibf_graph.root.favourite_child_idx.value(), 0);
    ASSERT_EQ(hibf_graph.root.remaining_records.size(), 2);
    EXPECT_EQ(hibf_graph.root.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 2, 1, 1}));
    EXPECT_EQ(hibf_graph.root.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{}, 3, 2, 0}));

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
    EXPECT_EQ(root_child0.favourite_child_idx.value(), 1);
    ASSERT_EQ(root_child0.remaining_records.size(), 3);
    EXPECT_EQ(root_child0.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0}, 2, 1, 7}));
    EXPECT_EQ(root_child0.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0}, 3, 1, 6}));
    EXPECT_EQ(root_child0.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0}, 4, 1, 5}));

    // child 1 of root
    EXPECT_EQ(root_child1.parent_bin_index, 1);
    EXPECT_EQ(root_child1.max_bin_index, 26);
    EXPECT_EQ(root_child1.number_of_technical_bins, 64);
    EXPECT_EQ(root_child1.favourite_child_idx.has_value(), false);
    ASSERT_EQ(root_child1.remaining_records.size(), 3);
    EXPECT_EQ(root_child1.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 26, 19, 3}));
    EXPECT_EQ(root_child1.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 0, 26, 2}));
    EXPECT_EQ(root_child1.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{1}, 45, 19, 4}));

    // child 0 of child 0
    EXPECT_EQ(c0_child0.parent_bin_index, 0);
    EXPECT_EQ(c0_child0.max_bin_index, 4);
    EXPECT_EQ(c0_child0.number_of_technical_bins, 5);
    EXPECT_EQ(c0_child0.favourite_child_idx.has_value(), false);
    ASSERT_EQ(c0_child0.remaining_records.size(), 4);
    EXPECT_EQ(c0_child0.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0}, 4, 1, 11}));
    EXPECT_EQ(c0_child0.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0}, 1, 1, 14}));
    EXPECT_EQ(c0_child0.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0}, 2, 1, 13}));
    EXPECT_EQ(c0_child0.remaining_records[3],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0}, 3, 1, 12}));

    // child 1 of child 0
    EXPECT_EQ(c0_child1.parent_bin_index, 1);
    EXPECT_EQ(c0_child1.max_bin_index, 34);
    EXPECT_EQ(c0_child1.number_of_technical_bins, 64);
    EXPECT_EQ(c0_child1.favourite_child_idx.has_value(), false);
    ASSERT_EQ(c0_child1.remaining_records.size(), 3);
    EXPECT_EQ(c0_child1.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 1}, 34, 15, 9}));
    EXPECT_EQ(c0_child1.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 1}, 0, 34, 8}));
    EXPECT_EQ(c0_child1.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 1}, 49, 15, 10}));

    // child 0 of child 0 of child 0
    EXPECT_EQ(c0_c0_child0.parent_bin_index, 0);
    EXPECT_EQ(c0_c0_child0.max_bin_index, 30);
    EXPECT_EQ(c0_c0_child0.number_of_technical_bins, 64);
    EXPECT_EQ(c0_c0_child0.favourite_child_idx.has_value(), false);
    ASSERT_EQ(c0_c0_child0.remaining_records.size(), 5);
    EXPECT_EQ(c0_c0_child0.remaining_records[0],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0, 0}, 30, 11, 16}));
    EXPECT_EQ(c0_c0_child0.remaining_records[1],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0, 0}, 0, 30, 15}));
    EXPECT_EQ(c0_c0_child0.remaining_records[2],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0, 0}, 41, 11, 17}));
    EXPECT_EQ(c0_c0_child0.remaining_records[3],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0, 0}, 52, 6, 18}));
    EXPECT_EQ(c0_c0_child0.remaining_records[4],
              (seqan::hibf::layout::layout::user_bin{std::vector<size_t>{0, 0, 0}, 58, 6, 19}));
}
