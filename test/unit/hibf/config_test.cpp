// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, TestPartResult, AssertionResult, Test, EXPECT_EQ, Capture...

#include <cstddef>    // for size_t
#include <functional> // for function
#include <sstream>    // for basic_stringstream, stringstream
#include <stdexcept>  // for invalid_argument
#include <string>     // for allocator, char_traits, string
#include <utility>    // for move

#include <hibf/config.hpp>                   // for config, insert_iterator
#include <hibf/misc/add_empty_bins.hpp>      // for add_empty_bins
#include <hibf/misc/subtract_empty_bins.hpp> // for subtract_empty_bins
#include <hibf/test/cereal.hpp>              // for test_serialisation
#include <hibf/test/expect_throw_msg.hpp>    // for EXPECT_THROW_MSG

TEST(config_test, write_to)
{
    std::stringstream ss{};

    seqan::hibf::config configuration;

    configuration.number_of_user_bins = 123456789;
    configuration.number_of_hash_functions = 4;
    configuration.maximum_fpr = 0.0001;
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
                                    "@        \"version\": 2,\n"
                                    "@        \"number_of_user_bins\": 123456789,\n"
                                    "@        \"number_of_hash_functions\": 4,\n"
                                    "@        \"maximum_fpr\": 0.0001,\n"
                                    "@        \"relaxed_fpr\": 0.3,\n"
                                    "@        \"threads\": 31,\n"
                                    "@        \"sketch_bits\": 8,\n"
                                    "@        \"tmax\": 128,\n"
                                    "@        \"empty_bin_fraction\": 0.0,\n"
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
                         "@        \"version\": 2,\n"
                         "@        \"number_of_user_bins\": 123456789,\n"
                         "@        \"number_of_hash_functions\": 4,\n"
                         "@        \"maximum_fpr\": 0.0001,\n"
                         "@        \"relaxed_fpr\": 0.3,\n"
                         "@        \"threads\": 31,\n"
                         "@        \"sketch_bits\": 8,\n"
                         "@        \"tmax\": 128,\n"
                         "@        \"empty_bin_fraction\": 0.5,\n"
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
    EXPECT_EQ(configuration.maximum_fpr, 0.0001);
    EXPECT_EQ(configuration.relaxed_fpr, 0.3);
    EXPECT_EQ(configuration.threads, 31);
    EXPECT_EQ(configuration.sketch_bits, 8);
    EXPECT_EQ(configuration.tmax, 128);
    EXPECT_EQ(configuration.empty_bin_fraction, 0.5);
    EXPECT_EQ(configuration.alpha, 1.0);
    EXPECT_EQ(configuration.max_rearrangement_ratio, 0.333);
    EXPECT_EQ(configuration.disable_estimate_union, true);
    EXPECT_EQ(configuration.disable_rearrangement, false);
}

TEST(config_test, read_from_v1)
{
    std::stringstream ss{"@HIBF_CONFIG\n"
                         "@{\n"
                         "@    \"hibf_config\": {\n"
                         "@        \"version\": 1,\n"
                         "@        \"number_of_user_bins\": 123456789,\n"
                         "@        \"number_of_hash_functions\": 4,\n"
                         "@        \"maximum_fpr\": 0.0001,\n"
                         "@        \"relaxed_fpr\": 0.3,\n"
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
    EXPECT_EQ(configuration.maximum_fpr, 0.0001);
    EXPECT_EQ(configuration.relaxed_fpr, 0.3);
    EXPECT_EQ(configuration.threads, 31);
    EXPECT_EQ(configuration.sketch_bits, 8);
    EXPECT_EQ(configuration.tmax, 128);
    EXPECT_EQ(configuration.empty_bin_fraction, 0.0);
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
                         "@        \"maximum_fpr\": 0.0001,\n"
                         "@        \"relaxed_fpr\": 0.3,\n"
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
    EXPECT_EQ(configuration.maximum_fpr, 0.0001);
    EXPECT_EQ(configuration.relaxed_fpr, 0.3);
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
    auto dummy_input_fn = [](size_t const, seqan::hibf::insert_iterator) {};

    // input_fn is not set
    {
        seqan::hibf::config configuration{};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] You did not set the required config::input_fn.");
    }

    // number_of_user_bins cannot be 0, bin_kind::merged (18'446'744'073'709'551'615ULL),
    // or bin_kind::deleted (bin_kind::merged - 1)
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] You did not set the required config::number_of_user_bins.");

        configuration.number_of_user_bins = 18'446'744'073'709'551'615ULL;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] The maximum possible config::number_of_user_bins "
                         "is 18446744073709551613.");

        configuration.number_of_user_bins = 18'446'744'073'709'551'614ULL;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] The maximum possible config::number_of_user_bins "
                         "is 18446744073709551613.");
    }

    // number_of_hash_functions must be in [1,5]
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .number_of_hash_functions = 0u};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::number_of_hash_functions must be in [1,5].");

        configuration.number_of_hash_functions = 6u;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::number_of_hash_functions must be in [1,5].");
    }

    // maximum_fpr must be in (0.0,1.0)
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 1u, .maximum_fpr = 0.0};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::maximum_fpr must be in (0.0,1.0).");

        configuration.maximum_fpr = 1.0;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::maximum_fpr must be in (0.0,1.0).");
    }

    // relaxed_fpr must be in (0.0,1.0)
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 1u, .relaxed_fpr = 0.0};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::relaxed_fpr must be in (0.0,1.0).");

        configuration.relaxed_fpr = 1.0;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::relaxed_fpr must be in (0.0,1.0).");
    }

    // relaxed_fpr must equal to or greater than maximum_fpr
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .maximum_fpr = 0.3,
                                          .relaxed_fpr = 0.2};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::relaxed_fpr must be "
                         "greater than or equal to config::maximum_fpr.");
    }

    // threads cannot be 0
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 1u, .threads = 0u};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::threads must be greater than 0.");
    }

    // sketch_bits must be in [5,32]
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 1u, .sketch_bits = 4u};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::sketch_bits must be in [5,32].");

        configuration.sketch_bits = 33u;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::sketch_bits must be in [5,32].");
    }

    // Set default tmax
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 4'286'582'784ULL};
        EXPECT_EQ(configuration.tmax, 0u);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_EQ(configuration.tmax, 65472u);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }

    // Given tmax OK
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .tmax = 18'446'744'073'709'551'552ULL};

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();

        EXPECT_EQ(configuration.number_of_user_bins, 1u);
        EXPECT_EQ(configuration.tmax, 18'446'744'073'709'551'552ULL);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }

    // Given tmax too big
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .tmax = 18'446'744'073'709'551'553ULL};

        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] The maximum possible config::tmax "
                         "is 18446744073709551552.");
    }

    // Given tmax is not a multiple of 64
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 32u, .tmax = 32u};
        EXPECT_EQ(configuration.tmax, 32u);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_EQ(configuration.tmax, 64u);
        EXPECT_EQ((testing::internal::GetCapturedStderr()),
                  "[HIBF CONFIG WARNING]: Your requested number of technical bins was not a multiple of 64. Due to the "
                  "architecture of the HIBF, it will use up space equal to the next multiple of 64 anyway, so we "
                  "increased your number of technical bins to 64.\n");
    }

    // empty_bin_fraction must be in [0.0,1.0)
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .empty_bin_fraction = -0.1};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::empty_bin_fraction must be in [0.0,1.0).");

        configuration.empty_bin_fraction = 1.0;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::empty_bin_fraction must be in [0.0,1.0).");
    }

    // alpha must be positive
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn, .number_of_user_bins = 1u, .alpha = -0.1};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::alpha must be positive.");
    }

    // max_rearrangement_ratio must be in [0.0,1.0]
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .max_rearrangement_ratio = -0.1};
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::max_rearrangement_ratio must be in [0.0,1.0].");

        configuration.max_rearrangement_ratio = 1.1;
        EXPECT_THROW_MSG(configuration.validate_and_set_defaults(),
                         std::invalid_argument,
                         "[HIBF CONFIG ERROR] config::max_rearrangement_ratio must be in [0.0,1.0].");
    }

    // Set disable_rearrangement if disable_estimate_union is set
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .disable_estimate_union = true};
        EXPECT_TRUE(configuration.disable_estimate_union);
        EXPECT_FALSE(configuration.disable_rearrangement);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_TRUE(configuration.disable_estimate_union);
        EXPECT_TRUE(configuration.disable_rearrangement);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }

    // Set disable_rearrangement if max_rearrangement_ratio is 0.0
    {
        seqan::hibf::config configuration{.input_fn = dummy_input_fn,
                                          .number_of_user_bins = 1u,
                                          .max_rearrangement_ratio = 0.0};
        EXPECT_EQ(configuration.max_rearrangement_ratio, 0.0);
        EXPECT_FALSE(configuration.disable_estimate_union);
        EXPECT_FALSE(configuration.disable_rearrangement);

        testing::internal::CaptureStderr();
        configuration.validate_and_set_defaults();
        EXPECT_EQ(configuration.max_rearrangement_ratio, 0.0);
        EXPECT_FALSE(configuration.disable_estimate_union);
        EXPECT_TRUE(configuration.disable_rearrangement);
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }
}

TEST(config_test, empty_bin_fraction)
{
    seqan::hibf::config config{.input_fn = [](size_t const, seqan::hibf::insert_iterator) {},
                               .number_of_user_bins = 1u};

    for (size_t top_level_tmax = 64; top_level_tmax <= 1024; top_level_tmax += 64)
    {
        config.tmax = top_level_tmax;

        // Some deviation for high fractions
        for (double fraction = 0.00; fraction <= 0.97; fraction += 0.01)
        {
            size_t const top_level_subtracted = seqan::hibf::subtract_empty_bins(top_level_tmax, fraction);
            double const adjusted_fraction = 1.0 - static_cast<double>(top_level_subtracted) / top_level_tmax;

            config.empty_bin_fraction = fraction;
            config.validate_and_set_defaults();
            EXPECT_EQ(config.empty_bin_fraction, adjusted_fraction);

            // Roundtrip test
            for (size_t lower_level_tmax = 64; lower_level_tmax <= top_level_tmax; lower_level_tmax += 64)
            {
                size_t const lower_level_subtracted =
                    seqan::hibf::subtract_empty_bins(lower_level_tmax, adjusted_fraction);
                size_t const lower_level_added = seqan::hibf::add_empty_bins(lower_level_subtracted, adjusted_fraction);

                ASSERT_EQ(seqan::hibf::next_multiple_of_64(lower_level_added), lower_level_tmax)
                    << "top_level_tmax: " << top_level_tmax << '\n'
                    << "lower_level_tmax: " << lower_level_tmax << '\n'
                    << "fraction: " << fraction << '\n'
                    << "top_level_subtracted: " << top_level_subtracted << '\n'
                    << "adjusted_fraction: " << adjusted_fraction << '\n'
                    << "lower_level_subtracted: " << lower_level_subtracted << '\n'
                    << "lower_level_added: " << lower_level_added << '\n';
            }
        }
    }
}

TEST(config_test, serialisation)
{
    seqan::hibf::config config{.number_of_user_bins = 123456789,
                               .number_of_hash_functions = 4,
                               .maximum_fpr = 0.0001,
                               .threads = 31,
                               .sketch_bits = 8,
                               .tmax = 128,
                               .alpha = 1.0,
                               .max_rearrangement_ratio = 0.333,
                               .disable_estimate_union = true,
                               .disable_rearrangement = false};
    seqan::hibf::test::test_serialisation(std::move(config));
}
