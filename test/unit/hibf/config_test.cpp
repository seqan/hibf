#include <gtest/gtest.h> // for Message, TestPartResult, Test, EXPECT_EQ, AssertionResult, CaptureStderr, GetCapt...

#include <sstream>   // for basic_stringstream, stringstream
#include <stdexcept> // for invalid_argument
#include <string>    // for allocator, string

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
                                    "@        \"disable_rearrangement\": false\n"
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
                         "@        \"disable_rearrangement\": false\n"
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
                         "@        \"disable_rearrangement\": false\n"
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
}

TEST(config_test, validate_and_set_defaults)
{
    // Set disable_rearrangement if disable_estimate_union is set
    {
        seqan::hibf::config configuration{.number_of_user_bins = 1u, .disable_estimate_union = true};
        EXPECT_TRUE(configuration.disable_estimate_union);
        EXPECT_FALSE(configuration.disable_rearrangement);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_TRUE(configuration.disable_estimate_union);
        EXPECT_TRUE(configuration.disable_rearrangement);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }

    // number_of_user_bins cannot be 0
    {
        seqan::hibf::config configuration{};

        try
        {
            configuration.validate_and_set_defaults();
            FAIL();
        }
        catch (std::invalid_argument const & exception)
        {
            EXPECT_STREQ("[HIBF CONFIG ERROR] You did not set the required config::number_of_user_bins.",
                         exception.what());
        }
        catch (...)
        {
            FAIL();
        }
    }

    // Set default tmax
    {
        seqan::hibf::config configuration{.number_of_user_bins = 4'286'582'784ULL};
        EXPECT_EQ(configuration.tmax, 0u);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_EQ(configuration.tmax, 65472u);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }

    // Cannot set default tmax
    {
        seqan::hibf::config configuration{.number_of_user_bins = 4'286'582'785ULL};
        EXPECT_EQ(configuration.tmax, 0u);

        try
        {
            configuration.validate_and_set_defaults();
            FAIL();
        }
        catch (std::invalid_argument const & exception)
        {
            EXPECT_STREQ(
                "[HIBF CONFIG ERROR] Too many user bins to compute a default tmax. Please set a tmax manually.",
                exception.what());
        }
        catch (...)
        {
            FAIL();
        }
    }

    // Given tmax OK
    {
        seqan::hibf::config configuration{.number_of_user_bins = 1u, .tmax = 65472u};

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();

        EXPECT_EQ(configuration.number_of_user_bins, 1u);
        EXPECT_EQ(configuration.tmax, 65472u);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }

    // Given tmax too big
    {
        seqan::hibf::config configuration{.number_of_user_bins = 1u, .tmax = 65473u};

        try
        {
            configuration.validate_and_set_defaults();
            FAIL();
        }
        catch (std::invalid_argument const & exception)
        {
            EXPECT_STREQ("[HIBF CONFIG ERROR] The maximum possible tmax is 65472.", exception.what());
        }
        catch (...)
        {
            FAIL();
        }
    }

    // Given tmax is not a multiple of 64
    {
        seqan::hibf::config configuration{.number_of_user_bins = 32u, .tmax = 32u};
        EXPECT_EQ(configuration.tmax, 32u);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_EQ(configuration.tmax, 64u);
        EXPECT_EQ((testing::internal::GetCapturedStderr()),
                  "[HIBF CONFIG WARNING]: Your requested number of technical bins was not a multiple of 64. Due to the "
                  "architecture of the HIBF, it will use up space equal to the next multiple of 64 anyway, so we "
                  "increased your number of technical bins to 64.\n");
    }
}
