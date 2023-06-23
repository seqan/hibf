#include <gtest/gtest.h>

#include <hibf/detail/build/parse_chopper_pack_header.hpp>

TEST(parse_chopper_pack_header_test, basic)
{
    std::stringstream chopper_pack_file
    {
     R"expected_layout(##CONFIG:
##{
##    "config": {
##        "version": 2,
##        "data_file": {
##            "value0": ""
##        },
##        "debug": false,
##        "sketch_directory": {
##            "value0": ""
##        },
##        "k": 19,
##        "sketch_bits": 12,
##        "disable_sketch_output": true,
##        "precomputed_files": false,
##        "output_filename": {
##            "value0": "dummy.layout"
##        },
##        "tmax": 128,
##        "num_hash_functions": 2,
##        "false_positive_rate": 0.05,
##        "alpha": 1.2,
##        "max_rearrangement_ratio": 0.5,
##        "threads": 1,
##        "disable_estimate_union": true,
##        "disable_rearrangement": true,
##        "determine_best_tmax": false,
##        "force_all_binnings": false
##    }
##}
##ENDCONFIG
#HIGH_LEVEL_IBF max_bin_id:14
#MERGED_BIN_14 max_bin_id:0
#MERGED_BIN_15 max_bin_id:0
#FILES	BIN_INDICES	NUMBER_OF_BINS
seq0	0	1
seq19	1	1
seq18	2	1
seq17	3	1
seq95	126	2
)expected_layout"
    };
    hibf::configuration chopper_config;
    hibf::layout::layout hibf_layout;

    hibf::parse_chopper_pack_header(chopper_pack_file, chopper_config, hibf_layout);

    // check chopper_configuration contents:
    EXPECT_EQ(chopper_config.debug, false);
    EXPECT_TRUE(chopper_config.sketch_directory.empty());
    EXPECT_EQ(chopper_config.k, 19);
    EXPECT_EQ(chopper_config.sketch_bits, 12);
    EXPECT_EQ(chopper_config.disable_sketch_output, true);
    EXPECT_EQ(chopper_config.precomputed_files, false);
    EXPECT_EQ(chopper_config.output_filename, "dummy.layout");
    EXPECT_EQ(chopper_config.tmax, 128);
    EXPECT_EQ(chopper_config.num_hash_functions, 2);
    EXPECT_EQ(chopper_config.false_positive_rate, 0.05);
    EXPECT_EQ(chopper_config.alpha, 1.2);
    EXPECT_EQ(chopper_config.max_rearrangement_ratio, 0.5);
    EXPECT_EQ(chopper_config.threads, 1u);
    EXPECT_EQ(chopper_config.disable_estimate_union, true);
    EXPECT_EQ(chopper_config.disable_rearrangement, true);
    EXPECT_EQ(chopper_config.determine_best_tmax, false);
    EXPECT_EQ(chopper_config.force_all_binnings, false);
    EXPECT_EQ(chopper_config.output_verbose_statistics, false);

    // check max bins
    ASSERT_EQ(hibf_layout.max_bins.size(), 2u);
    EXPECT_EQ(hibf_layout.max_bins[0], (hibf::layout::layout::max_bin{{14}, 0}));
    EXPECT_EQ(hibf_layout.max_bins[1], (hibf::layout::layout::max_bin{{15}, 0}));
    EXPECT_EQ(hibf_layout.top_level_max_bin_id, 14);
}
