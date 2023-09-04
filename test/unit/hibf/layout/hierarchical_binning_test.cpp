#include <gtest/gtest.h> // for Test, Message, TestPartResult, TestInfo, EXPECT_EQ

#include <cstddef> // for size_t
#include <vector>  // for vector, allocator

#include <hibf/config.hpp>                        // for config
#include <hibf/layout/compute_fpr_correction.hpp> // for compute_fpr_correction
#include <hibf/layout/data_store.hpp>             // for data_store
#include <hibf/layout/hierarchical_binning.hpp>   // for hierarchical_binning
#include <hibf/layout/layout.hpp>                 // for layout
#include <hibf/test/expect_range_eq.hpp>          // for expect_range_eq, EXPECT_RANGE_EQ

TEST(hierarchical_binning_test, small_example)
{
    seqan::hibf::config config;
    config.tmax = 4;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{500, 1000, 500, 500, 500, 500, 500, 500};

    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});
    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:3

    std::vector<seqan::hibf::layout::layout::max_bin> expected_max_bins{{{1}, 22}, {{2}, 22}};

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{}, 0, 1, 7},
                                                                          {{1}, 0, 22, 4},
                                                                          {{1}, 22, 21, 5},
                                                                          {{1}, 43, 21, 6},
                                                                          {{2}, 0, 22, 0},
                                                                          {{2}, 22, 21, 2},
                                                                          {{2}, 43, 21, 3},
                                                                          {{}, 3, 1, 1}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, another_example)
{
    seqan::hibf::config config;
    config.tmax = 5;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{50, 1000, 1000, 50, 5, 10, 10, 5};
    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});

    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<seqan::hibf::layout::layout::max_bin> expected_max_bins{{{0, 0}, 56}, {{0}, 0}};

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{0, 0}, 0, 42, 6},
                                                                          {{0, 0}, 42, 14, 5},
                                                                          {{0, 0}, 56, 4, 7},
                                                                          {{0, 0}, 60, 4, 4},
                                                                          {{0}, 1, 2, 0},
                                                                          {{0}, 3, 2, 3},
                                                                          {{}, 1, 2, 2},
                                                                          {{}, 3, 2, 1}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, high_level_max_bin_id_is_0)
{
    seqan::hibf::config config;
    config.tmax = 4;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{500, 500, 500, 500};
    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});

    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{}, 0, 1, 3},
                                                                          {{}, 1, 1, 2},
                                                                          {{}, 2, 1, 1},
                                                                          {{}, 3, 1, 0}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, knuts_example)
{
    seqan::hibf::config config;
    config.alpha = 1;
    config.tmax = 5;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{60, 600, 1000, 800, 800};
    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});

    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u);

    std::vector<seqan::hibf::layout::layout::max_bin> expected_max_bins{{{0}, 63}};

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{0}, 0, 63, 1},
                                                                          {{0}, 63, 1, 0},
                                                                          {{}, 1, 1, 4},
                                                                          {{}, 2, 1, 3},
                                                                          {{}, 3, 2, 2}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, four_level_hibf)
{
    seqan::hibf::config config;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{11090, 5080, 3040, 1020, 510, 500};
    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});

    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1u); // #HIGH_LEVEL_IBF max_bin_id:1

    std::vector<seqan::hibf::layout::layout::max_bin> expected_max_bins{{{0, 0, 0, 0}, 33},
                                                                        {{0, 0, 0}, 1},
                                                                        {{0, 0}, 1},
                                                                        {{0}, 1}};

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{0, 0, 0, 0}, 0, 33, 4},
                                                                          {{0, 0, 0, 0}, 33, 31, 5},
                                                                          {{0, 0, 0}, 1, 1, 3},
                                                                          {{0, 0}, 1, 1, 2},
                                                                          {{0}, 1, 1, 1},
                                                                          {{}, 1, 1, 0}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin)
{
    seqan::hibf::config config;
    config.alpha = 1;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{500, 500, 500, 500};
    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});

    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::vector<seqan::hibf::layout::layout::max_bin> expected_max_bins{{{0}, 0}, {{1}, 0}};

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{0}, 0, 32, 2},
                                                                          {{0}, 32, 32, 3},
                                                                          {{1}, 0, 32, 0},
                                                                          {{1}, 32, 32, 1}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}

TEST(hierarchical_binning_test, tb0_is_a_merged_bin_and_leads_to_recursive_call)
{
    seqan::hibf::config config;
    config.alpha = 1;
    config.tmax = 2;
    config.disable_estimate_union = true; // also disables rearrangement

    seqan::hibf::layout::layout hibf_layout{};
    std::vector<size_t> kmer_counts{500, 500, 500, 500, 500, 500, 500, 500};
    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout, .kmer_counts = &kmer_counts};

    data.fpr_correction =
        seqan::hibf::layout::compute_fpr_correction({.fpr = 0.05, .hash_count = 2, .t_max = config.tmax});

    seqan::hibf::layout::hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 0u);

    std::vector<seqan::hibf::layout::layout::max_bin> expected_max_bins{{{0, 0}, 0},
                                                                        {{0, 1}, 0},
                                                                        {{0}, 0},
                                                                        {{1, 0}, 0},
                                                                        {{1, 1}, 0},
                                                                        {{1}, 0}};

    std::vector<seqan::hibf::layout::layout::user_bin> expected_user_bins{{{0, 0}, 0, 32, 5},
                                                                          {{0, 0}, 32, 32, 4},
                                                                          {{0, 1}, 0, 32, 7},
                                                                          {{0, 1}, 32, 32, 6},
                                                                          {{1, 0}, 0, 32, 1},
                                                                          {{1, 0}, 32, 32, 0},
                                                                          {{1, 1}, 0, 32, 3},
                                                                          {{1, 1}, 32, 32, 2}};

    EXPECT_RANGE_EQ(hibf_layout.max_bins, expected_max_bins);
    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected_user_bins);
}
