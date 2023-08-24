#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, Message, TEST, TestPartResult

#include <cstddef>     // for size_t
#include <sstream>     // for operator<<, char_traits, basic_ostream, basic_stringstream, strings...
#include <string>      // for allocator, string
#include <string_view> // for operator<<
#include <vector>      // for vector

#include <hibf/config.hpp> // for config

TEST(config_test, write_to)
{
    std::stringstream ss{};

    seqan::hibf::config configuration;

    configuration.number_of_user_bins = 123456789;
    configuration.number_of_hash_functions = 4;
    configuration.maximum_false_positive_rate = 0.0001;
    configuration.threads = 31;
    configuration.sketch_bits = 8;
    configuration.tmax = 128;
    configuration.alpha = 1.0;
    configuration.max_rearrangement_ratio = 0.333;
    configuration.disable_estimate_union = true;
    configuration.disable_rearrangement = false;
    configuration.disable_cutoffs = false;

    configuration.write_to(ss);

    std::string const expected_file{"@HIBF_CONFIG\n"
                                    "@{\n"
                                    "@    \"hibf_config\": {\n"
                                    "@        \"version\": 1,\n"
                                    "@        \"number_of_user_bins\": 123456789,\n"
                                    "@        \"number_of_hash_functions\": 4,\n"
                                    "@        \"maximum_false_positive_rate\": 0.0001,\n"
                                    "@        \"threads\": 31,\n"
                                    "@        \"sketch_bits\": 8,\n"
                                    "@        \"tmax\": 128,\n"
                                    "@        \"alpha\": 1.0,\n"
                                    "@        \"max_rearrangement_ratio\": 0.333,\n"
                                    "@        \"disable_estimate_union\": true,\n"
                                    "@        \"disable_rearrangement\": false,\n"
                                    "@        \"disable_cutoffs\": false\n"
                                    "@    }\n"
                                    "@}\n"
                                    "@HIBF_CONFIG_END\n"};

    EXPECT_EQ(ss.str(), expected_file);
}

TEST(config_test, read_from)
{
    std::stringstream ss{"@HIBF_CONFIG\n"
                         "@{\n"
                         "@    \"hibf_config\": {\n"
                         "@        \"version\": 1,\n"
                         "@        \"number_of_user_bins\": 123456789,\n"
                         "@        \"number_of_hash_functions\": 4,\n"
                         "@        \"maximum_false_positive_rate\": 0.0001,\n"
                         "@        \"threads\": 31,\n"
                         "@        \"sketch_bits\": 8,\n"
                         "@        \"tmax\": 128,\n"
                         "@        \"alpha\": 1.0,\n"
                         "@        \"max_rearrangement_ratio\": 0.333,\n"
                         "@        \"disable_estimate_union\": true,\n"
                         "@        \"disable_rearrangement\": false,\n"
                         "@        \"disable_cutoffs\": false\n"
                         "@    }\n"
                         "@}\n"
                         "@HIBF_CONFIG_END\n"};

    seqan::hibf::config configuration;
    configuration.read_from(ss);

    EXPECT_EQ(configuration.number_of_user_bins, 123456789);
    EXPECT_EQ(configuration.number_of_hash_functions, 4);
    EXPECT_EQ(configuration.maximum_false_positive_rate, 0.0001);
    EXPECT_EQ(configuration.threads, 31);
    EXPECT_EQ(configuration.sketch_bits, 8);
    EXPECT_EQ(configuration.tmax, 128);
    EXPECT_EQ(configuration.alpha, 1.0);
    EXPECT_EQ(configuration.max_rearrangement_ratio, 0.333);
    EXPECT_EQ(configuration.disable_estimate_union, true);
    EXPECT_EQ(configuration.disable_rearrangement, false);
    EXPECT_EQ(configuration.disable_cutoffs, false);
}

TEST(config_test, read_from_with_more_meta)
{
    std::stringstream ss{"@blah some chopper stuff\n"
                         "@blah some chopper stuff\n"
                         "@blah some chopper stuff\n"
                         "@blah some chopper stuff\n"
                         "@blah some chopper stuff\n"
                         "@HIBF_CONFIG\n"
                         "@{\n"
                         "@    \"hibf_config\": {\n"
                         "@        \"version\": 1,\n"
                         "@        \"number_of_user_bins\": 123456789,\n"
                         "@        \"number_of_hash_functions\": 4,\n"
                         "@        \"maximum_false_positive_rate\": 0.0001,\n"
                         "@        \"threads\": 31,\n"
                         "@        \"sketch_bits\": 8,\n"
                         "@        \"tmax\": 128,\n"
                         "@        \"alpha\": 1.0,\n"
                         "@        \"max_rearrangement_ratio\": 0.333,\n"
                         "@        \"disable_estimate_union\": true,\n"
                         "@        \"disable_rearrangement\": false,\n"
                         "@        \"disable_cutoffs\": false\n"
                         "@    }\n"
                         "@}\n"
                         "@HIBF_CONFIG_END\n"};

    seqan::hibf::config configuration;
    configuration.read_from(ss);

    EXPECT_EQ(configuration.number_of_user_bins, 123456789);
    EXPECT_EQ(configuration.number_of_hash_functions, 4);
    EXPECT_EQ(configuration.maximum_false_positive_rate, 0.0001);
    EXPECT_EQ(configuration.threads, 31);
    EXPECT_EQ(configuration.sketch_bits, 8);
    EXPECT_EQ(configuration.tmax, 128);
    EXPECT_EQ(configuration.alpha, 1.0);
    EXPECT_EQ(configuration.max_rearrangement_ratio, 0.333);
    EXPECT_EQ(configuration.disable_estimate_union, true);
    EXPECT_EQ(configuration.disable_rearrangement, false);
    EXPECT_EQ(configuration.disable_cutoffs, false);
}
