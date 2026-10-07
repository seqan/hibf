// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, Message, AssertionResult, TestPartResult, CmpHelperEQ, CmpHelperEQFa...

#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <format>      // for format
#include <functional>  // for function
#include <limits>      // for numeric_limits
#include <optional>    // for optional
#include <ranges>      // for transform
#include <sstream>     // for basic_stringstream, operator<<, stringstream, basic_ios, basic_iostream
#include <stdexcept>   // for invalid_argument
#include <string>      // for char_traits, allocator, basic_string, string
#include <string_view> // for basic_string_view
#include <utility>     // for pair, swap
#include <vector>      // for vector

#include <hibf/config.hpp>                      // for config, insert_iterator
#include <hibf/layout/compute_layout.hpp>       // for compute_layout
#include <hibf/layout/layout.hpp>               // for layout, operator<<
#include <hibf/sketch/compute_sketches.hpp>     // for compute_sketches
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog
#include <hibf/test/expect_range_eq.hpp>        // for EXPECT_RANGE_EQ
#include <hibf/test/expect_throw_msg.hpp>       // for EXPECT_THROW_MSG

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
using diagnostic_t = layout_t::diagnostic;
using code_t = diagnostic_t::code;
using severity_t = diagnostic_t::severity;

seqan::hibf::config const valid_config{.number_of_user_bins = 6u, .tmax = 64u};

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

// Returns the result of `layout.validate` and all reported diagnostics.
std::pair<bool, std::vector<diagnostic_t>> validate(layout_t const & layout,
                                                    seqan::hibf::config const & config = valid_config)
{
    std::vector<diagnostic_t> diagnostics{};
    bool const valid = layout.validate(config,
                                       [&diagnostics](diagnostic_t const & diagnostic)
                                       {
                                           diagnostics.push_back(diagnostic);
                                       });
    return {valid, diagnostics};
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

TEST(layout_test, validate)
{
    layout_t const layout = valid_layout();

    auto const [valid, diagnostics] = validate(layout);
    EXPECT_TRUE(valid);
    EXPECT_TRUE(diagnostics.empty());

    // The user bins below a merged bin do not need to be consecutive.
    // Order: UB 5, UB 0, UB 3 (IBF 2;2), UB 2, UB 1 (IBF 2), UB 4 (IBF 2;2).
    layout_t interleaved = layout;
    std::swap(interleaved.user_bins[2], interleaved.user_bins[4]);
    EXPECT_TRUE(validate(interleaved).second.empty());
}

TEST(layout_test, format_diagnostic)
{
    diagnostic_t const note{.level = severity_t::note, .message = "a"};
    diagnostic_t const warning{.level = severity_t::warning, .message = "b"};
    diagnostic_t const error{.level = severity_t::error, .message = "c"};

    EXPECT_EQ(std::format("{}", note), "[HIBF LAYOUT NOTE] a");
    EXPECT_EQ(std::format("{}", warning), "[HIBF LAYOUT WARNING] b");
    EXPECT_EQ(std::format("{}", error), "[HIBF LAYOUT ERROR] c");
    EXPECT_EQ(std::format("{:}", error), "[HIBF LAYOUT ERROR] c");
    EXPECT_THROW_MSG(static_cast<void>(std::vformat("{:>30}", std::make_format_args(error))),
                     std::format_error,
                     "seqan::hibf::layout::layout::diagnostic does not support format specifications.");

    std::ostringstream stream{};
    stream << warning;
    EXPECT_EQ(stream.str(), "[HIBF LAYOUT WARNING] b");

    EXPECT_THROW_MSG(static_cast<void>(std::format("{}", diagnostic_t{.level = static_cast<severity_t>(3)})),
                     std::format_error,
                     "Invalid seqan::hibf::layout::layout::diagnostic::severity.");
}

TEST(layout_test, throw_on_error)
{
    EXPECT_THROW_MSG(layout_t::throw_on_error({.level = severity_t::error, .message = "c"}),
                     std::invalid_argument,
                     "[HIBF LAYOUT ERROR] c");

    // Warnings are printed; notes only without NDEBUG.
    testing::internal::CaptureStderr();
    layout_t::throw_on_error({.level = severity_t::warning, .message = "b"});
    layout_t::throw_on_error({.level = severity_t::note, .message = "a"});
    std::string expected_stderr{"[HIBF LAYOUT WARNING] b\n"};
#ifndef NDEBUG
    expected_stderr += "[HIBF LAYOUT NOTE] a\n";
#endif
    EXPECT_EQ(testing::internal::GetCapturedStderr(), expected_stderr);
}

TEST(layout_test, validate_errors)
{
    struct defect
    {
        std::string_view description;
        code_t expected;
        std::function<void(layout_t &)> apply;
    };

    // Each defect is applied to valid_layout() and must be the only reported diagnostic.
    // clang-format off
    std::vector<defect> const defects{
        {"no user bins", code_t::empty_layout,
         [](layout_t & l) { l.clear(); }},
        {"zero technical bins", code_t::no_technical_bins,
         [](layout_t & l) { l.user_bins[0].number_of_technical_bins = 0u; }},
        {"technical bin range overflows", code_t::technical_bin_overflow,
         [](layout_t & l) { l.user_bins[0] = {{}, std::numeric_limits<size_t>::max() - 64u, 2u, 5u}; }},
        {"more than 2^64 - 64 technical bins", code_t::technical_bin_overflow, // next_multiple_of_64 overflows
         [](layout_t & l) { l.user_bins[0] = {{}, std::numeric_limits<size_t>::max() - 63u, 1u, 5u}; }},
        {"merged bin index overflows", code_t::technical_bin_overflow, // graph.cpp computes merged bin + 1
         [](layout_t & l) { l.user_bins[4].previous_TB_indices[1] = std::numeric_limits<size_t>::max() - 63u;
                            l.user_bins[5].previous_TB_indices[1] = std::numeric_limits<size_t>::max() - 63u;
                            l.max_bins[1].previous_TB_indices[1] = std::numeric_limits<size_t>::max() - 63u; }},
        {"user bin index out of range", code_t::user_bin_out_of_range,
         [](layout_t & l) { l.user_bins[0].idx = 6u; }},
        {"duplicate user bin index", code_t::duplicate_user_bin,
         [](layout_t & l) { l.user_bins[0].idx = 0u; }},
        {"user bins overlap", code_t::overlapping_technical_bins, // UB 0 occupies TBs 0-1
         [](layout_t & l) { l.user_bins[0].storage_TB_id = 1u; }},
        {"user bin in merged bin", code_t::overlapping_technical_bins,
         [](layout_t & l) { l.user_bins[0].storage_TB_id = 2u; }},
        {"user bin passes through user bin", code_t::overlapping_technical_bins, // TB 1 of IBF 2 holds UB 1
         [](layout_t & l) { l.user_bins[4].previous_TB_indices = l.user_bins[5].previous_TB_indices = {2, 1};
                            l.max_bins[1].previous_TB_indices = {2, 1}; }},
        {"max bin for top-level IBF", code_t::top_level_max_bin_entry,
         [](layout_t & l) { l.max_bins.insert(l.max_bins.begin(), layout_t::max_bin{{}, 2}); }},
        {"max bin for non-existing IBF", code_t::max_bin_without_ibf,
         [](layout_t & l) { l.max_bins[1].previous_TB_indices = {2, 3}; }},
        {"max bin duplicated", code_t::duplicate_max_bin,
         [](layout_t & l) { l.max_bins[1] = l.max_bins[0]; }},
        {"max bin refers to technical bin beyond the IBF", code_t::invalid_max_bin,
         [](layout_t & l) { l.max_bins[0].id = 64u; }},
        {"max bin missing", code_t::missing_max_bin,
         [](layout_t & l) { l.max_bins.pop_back(); }},
        {"top-level max bin refers to technical bin beyond the IBF", code_t::invalid_max_bin,
         [](layout_t & l) { l.top_level_max_bin_id = 64u; }},
        {"top-level max bin refers to second technical bin of split bin", code_t::invalid_max_bin,
         [](layout_t & l) { l.top_level_max_bin_id = 1u; }}};
    // clang-format on

    for (auto const & [description, expected, apply] : defects)
    {
        SCOPED_TRACE(description);
        layout_t layout = valid_layout();
        apply(layout);

        auto const [valid, diagnostics] = validate(layout);
        EXPECT_FALSE(valid);
        EXPECT_EQ(diagnostics.size(), 1u);
        if (diagnostics.empty())
            continue; // not ASSERT_EQ: it would skip the remaining defects
        EXPECT_EQ(diagnostics[0].level, severity_t::error);
        EXPECT_EQ(diagnostics[0].what, expected);

        EXPECT_FALSE(layout.validate(valid_config, {}));
        EXPECT_THROW_MSG(layout.validate(valid_config),
                         std::invalid_argument,
                         "[HIBF LAYOUT ERROR] " + diagnostics[0].message);
    }
}

TEST(layout_test, validate_error_details)
{
    // The first user bin index that is not present is reported.
    layout_t layout = valid_layout();
    layout.user_bins[4].idx = 6u; // UB 3 -> UB 6
    {
        auto const [valid, diagnostics] = validate(layout, seqan::hibf::config{.number_of_user_bins = 7u});
        ASSERT_EQ(diagnostics.size(), 1u);
        EXPECT_EQ(diagnostics[0].what, code_t::missing_user_bin);
        EXPECT_EQ(diagnostics[0].user_bin, 3u);
        EXPECT_FALSE(diagnostics[0].ibf.has_value()); // concerns no IBF
    }
    EXPECT_THROW_MSG(layout.validate(seqan::hibf::config{.number_of_user_bins = 7u}),
                     std::invalid_argument,
                     "[HIBF LAYOUT ERROR] User bin 3 is missing from the layout.");

    // The IBF, the user bin and the technical bin of an overlap are reported.
    layout = valid_layout();
    layout.user_bins[3].storage_TB_id = 0u; // UB 2 -> TB 0 of IBF 2, which holds UB 1
    {
        auto const [valid, diagnostics] = validate(layout);
        ASSERT_EQ(diagnostics.size(), 1u);
        EXPECT_EQ(diagnostics[0].ibf, std::vector<size_t>{2u});
        EXPECT_EQ(diagnostics[0].user_bin, 2u);
        EXPECT_EQ(diagnostics[0].technical_bin, 0u);
        EXPECT_EQ(diagnostics[0].message, "In IBF 2, technical bin 0 is used by user bin 1 and by user bin 2.");
    }

    // A missing max bin names the IBF.
    layout = valid_layout();
    layout.max_bins.pop_back();
    {
        auto const [valid, diagnostics] = validate(layout);
        ASSERT_EQ(diagnostics.size(), 1u);
        EXPECT_EQ(diagnostics[0].ibf, (std::vector<size_t>{2u, 2u}));
        EXPECT_EQ(diagnostics[0].message, "The max bins (\"#LOWER_LEVEL_IBF\") contain no entry for IBF 2;2.");
    }

    // An overflowing merged bin is reported as technical bin of the IBF that contains it.
    size_t const overflowing = std::numeric_limits<size_t>::max() - 63u; // IBF 2;<overflowing> has 2^64 - 63 TBs
    layout = valid_layout();
    layout.user_bins[4].previous_TB_indices[1] = layout.user_bins[5].previous_TB_indices[1] = overflowing;
    layout.max_bins[1].previous_TB_indices[1] = overflowing;
    {
        auto const [valid, diagnostics] = validate(layout);
        ASSERT_EQ(diagnostics.size(), 1u);
        EXPECT_EQ(diagnostics[0].ibf, std::vector<size_t>{2u});
        EXPECT_EQ(diagnostics[0].user_bin, 3u);
        EXPECT_EQ(diagnostics[0].technical_bin, overflowing);
        EXPECT_EQ(diagnostics[0].message,
                  "In IBF 2, the merged bin on the path of user bin 3 exceeds the maximum number of technical bins "
                  "of an IBF.");
    }

    // 2^64 - 64 technical bins are supported: next_multiple_of_64 does not overflow.
    layout = valid_layout();
    layout.user_bins[0] = {{}, overflowing - 1u, 1u, 5u}; // UB 5 in the last supported technical bin
    EXPECT_TRUE(layout.validate(valid_config, {}));

    layout = valid_layout();
    layout.max_bins.emplace_back(std::vector<size_t>{}, 2u);
    EXPECT_THROW_MSG(layout.validate(valid_config),
                     std::invalid_argument,
                     "[HIBF LAYOUT ERROR] The max bins contain an entry for the Root-IBF. The Root-IBF's max bin is "
                     "given by top_level_max_bin_id (\"#TOP_LEVEL_IBF\"), not by the max bins "
                     "(\"#LOWER_LEVEL_IBF\").");
}

TEST(layout_test, validate_warnings_and_notes)
{
    // Technical bins beyond tmax: UB 5 moves to TBs 64-127. The Root-IBF then also has empty technical bins 3 to 63.
    layout_t layout = valid_layout();
    layout.user_bins[0] = {{}, 64u, 64u, 5u};
    {
        auto const [valid, diagnostics] = validate(layout);
        EXPECT_TRUE(valid);
        ASSERT_EQ(diagnostics.size(), 2u);
        EXPECT_EQ(diagnostics[0].level, severity_t::warning);
        EXPECT_EQ(diagnostics[0].what, code_t::technical_bin_exceeds_tmax);
        EXPECT_EQ(diagnostics[0].ibf, std::vector<size_t>{});
        EXPECT_EQ(diagnostics[0].technical_bin, 127u);
        EXPECT_EQ(diagnostics[1].level, severity_t::note);
        EXPECT_EQ(diagnostics[1].what, code_t::empty_technical_bins);
        EXPECT_EQ(diagnostics[1].technical_bin, 3u);
        EXPECT_EQ(diagnostics[1].message,
                  "The Root-IBF uses technical bins 0-127, but 61 of them are empty; the first one is 3.");
    }

    // tmax is not checked if it is 0.
    {
        auto const [valid, diagnostics] = validate(layout, seqan::hibf::config{.number_of_user_bins = 6u});
        EXPECT_TRUE(valid);
        ASSERT_EQ(diagnostics.size(), 1u);
        EXPECT_EQ(diagnostics[0].what, code_t::empty_technical_bins);
    }

    // tmax is rounded up to a multiple of 64. The message names a tmax that is not a multiple of 64.
    layout.user_bins[0] = {{}, 128u, 64u, 5u}; // UB 5 in TBs 128-191
    {
        auto const [valid, diagnostics] =
            validate(layout, seqan::hibf::config{.number_of_user_bins = 6u, .tmax = 100u});
        ASSERT_EQ(diagnostics.size(), 2u);
        EXPECT_EQ(diagnostics[0].message,
                  "The Root-IBF uses technical bin 191, but tmax 100 (rounded up to a multiple of 64) allows only 128 "
                  "technical bins.");
    }

    // A lower-level IBF with a single user bin: UB 4 moves from IBF 2;2 to TB 63 of IBF 2.
    layout = valid_layout();
    layout.user_bins[3].number_of_technical_bins = 60u; // UB 2 in TBs 3-62 of IBF 2
    layout.user_bins[5] = {{2u}, 63u, 1u, 4u};
    layout.user_bins[4].number_of_technical_bins = 64u; // UB 3 in TBs 0-63 of IBF 2;2
    layout.max_bins[1].id = 0u;
    {
        auto const [valid, diagnostics] = validate(layout);
        EXPECT_TRUE(valid);
        ASSERT_EQ(diagnostics.size(), 1u);
        EXPECT_EQ(diagnostics[0].level, severity_t::warning);
        EXPECT_EQ(diagnostics[0].what, code_t::single_bin_ibf);
        EXPECT_EQ(diagnostics[0].ibf, (std::vector<size_t>{2u, 2u}));
        EXPECT_EQ(diagnostics[0].user_bin, 3u);
    }

    // A lower-level IBF with a single merged bin: IBF 2 only contains merged bin 0; IBF 2;0 holds UBs 1-4 in TBs 0-63.
    // IBF 2 then also ends with 63 empty technical bins.
    layout = valid_layout();
    layout.max_bins[1] = {{2u, 0u}, 0u};
    for (size_t i = 2u; i < 6u; ++i)
        layout.user_bins[i] = {{2u, 0u}, i - 2u, (i == 5u) ? 61u : 1u, i - 1u};
    {
        auto const [valid, diagnostics] = validate(layout);
        EXPECT_TRUE(valid);
        ASSERT_EQ(diagnostics.size(), 2u);
        EXPECT_EQ(diagnostics[0].what, code_t::single_bin_ibf);
        EXPECT_EQ(diagnostics[0].ibf, std::vector<size_t>{2u});
        EXPECT_FALSE(diagnostics[0].user_bin.has_value());
        EXPECT_EQ(diagnostics[0].message,
                  "IBF 2 only contains merged bin 0. Moving the IBF below it up would save a level.");
        EXPECT_EQ(diagnostics[1].what, code_t::unexpected_empty_bins);
    }
}

// The technical bins of an IBF should be its used technical bins, followed by as many empty technical bins as
// empty_bin_fraction implies, i.e., none if it is 0. The build rounds the technical bins up to a multiple of 64.
TEST(layout_test, validate_empty_bins)
{
    // A Root-IBF with one user bin in every `stride`-th technical bin.
    auto root_ibf = [](size_t const user_bins, size_t const stride)
    {
        layout_t layout{};
        for (size_t i = 0u; i < user_bins; ++i)
            layout.user_bins.emplace_back(std::vector<size_t>{}, i * stride, 1u, i);
        return layout;
    };

    struct test_case
    {
        std::string_view description;
        size_t user_bins;
        size_t stride;
        double empty_bin_fraction;
        std::vector<code_t> expected;
    };

    // config::validate_and_set_defaults() sets empty_bin_fraction 0.1 to 1 - 57 / 64 = 0.109375 for tmax 64: 7 of 64
    // technical bins are empty. For 0.015625 = 1 / 64, 1 of 64 technical bins is empty.
    // clang-format off
    std::vector<test_case> const test_cases{
        {"empty bins without empty_bin_fraction", 57u, 1u, 0.0, {code_t::unexpected_empty_bins}},
        {"as many empty bins as expected", 57u, 1u, 0.109375, {}},
        {"one empty bin too many", 56u, 1u, 0.109375, {code_t::unexpected_empty_bins}},
        {"alternating", 32u, 2u, 0.015625, {code_t::empty_technical_bins}}, // TB 63 is the only trailing empty bin
        {"technical bins with empty bins do not fit into size_t", 2u, size_t{1} << 62, 0.6,
         {code_t::technical_bin_exceeds_tmax, code_t::empty_technical_bins}}};
    // clang-format on

    for (auto const & [description, user_bins, stride, empty_bin_fraction, expected] : test_cases)
    {
        SCOPED_TRACE(description);
        auto const [valid, diagnostics] =
            validate(root_ibf(user_bins, stride),
                     {.number_of_user_bins = user_bins, .tmax = 64u, .empty_bin_fraction = empty_bin_fraction});
        EXPECT_TRUE(valid);
        EXPECT_RANGE_EQ(diagnostics | std::views::transform(&diagnostic_t::what), expected);
    }

    // Alternating without empty_bin_fraction: Intermittent empty bins do not count as empty bins at the end.
    {
        auto const [valid, diagnostics] = validate(root_ibf(32u, 2u), {.number_of_user_bins = 32u, .tmax = 64u});
        ASSERT_EQ(diagnostics.size(), 2u);
        EXPECT_EQ(diagnostics[0].technical_bin, 62u);
        EXPECT_EQ(
            diagnostics[0].message,
            "The Root-IBF ends with 1 empty technical bin, but 0 are expected. There is a total of 64 technical bins "
            "and the empty_bin_fraction is 0.");
        EXPECT_EQ(diagnostics[1].technical_bin, 1u);
        EXPECT_EQ(diagnostics[1].message,
                  "The Root-IBF uses technical bins 0-62, but 31 of them are empty; the first one is 1.");
    }
    // 58 / (1 - 0.109375) > 64: The build uses 128 technical bins.
    {
        auto const [valid, diagnostics] =
            validate(root_ibf(58u, 1u), {.number_of_user_bins = 58u, .tmax = 64u, .empty_bin_fraction = 0.109375});
        ASSERT_EQ(diagnostics.size(), 2u);
        EXPECT_EQ(diagnostics[0].message,
                  "The Root-IBF uses technical bin 57. With the empty bins added by the build (empty_bin_fraction "
                  "0.109375), it exceeds the 64 technical bins tmax allows.");
        EXPECT_EQ(
            diagnostics[1].message,
            "The Root-IBF ends with 70 empty technical bins, but 14 are expected. There is a total of 128 technical "
            "bins and the empty_bin_fraction is 0.109375.");
    }
}

// The build sizes fpr_correction by the number of technical bins of the Root-IBF (next_multiple_of_64(#TBs) + 1 values)
// and indexes it with the number of technical bins of each IBF's max bin.
TEST(layout_test, validate_max_bin_split_beyond_root)
{
    // Root-IBF: UB 0 in TB 0, merged bin in TB 1 (max bin). fpr_correction has next_multiple_of_64(2) + 1 = 65 values.
    // IBF 1: UB 1 split into TBs 0-199 (max bin), UB 2 in TB 200. The build reads fpr_correction[200].
    layout_t layout{};
    layout.top_level_max_bin_id = 1u;
    layout.max_bins.emplace_back(std::vector<size_t>{1u}, 0u);
    layout.user_bins.emplace_back(std::vector<size_t>{}, 0u, 1u, 0u);
    layout.user_bins.emplace_back(std::vector<size_t>{1u}, 0u, 200u, 1u);
    layout.user_bins.emplace_back(std::vector<size_t>{1u}, 200u, 1u, 2u);

    seqan::hibf::config const config{.number_of_user_bins = 3u, .tmax = 256u};
    auto const [valid, diagnostics] = validate(layout, config);
    EXPECT_FALSE(valid);
    ASSERT_EQ(diagnostics.size(), 1u);
    EXPECT_EQ(diagnostics[0].what, code_t::max_bin_exceeds_root);
    EXPECT_EQ(diagnostics[0].ibf, std::vector<size_t>{1u});
    EXPECT_EQ(diagnostics[0].user_bin, 1u);
    EXPECT_EQ(diagnostics[0].message,
              "The max bin (\"fullest_technical_bin_idx:\") of IBF 1 spans 200 technical bins, but at most 64 "
              "are supported (the Root-IBF's 2 technical bins, rounded up to a multiple of 64).");
    EXPECT_THROW_MSG(layout.validate(config), std::invalid_argument, "[HIBF LAYOUT ERROR] " + diagnostics[0].message);
}

// Why using "IBF 0" for the top-level IBF would be ambiguous.
TEST(layout_test, validate_ibf_names)
{
    // Root-IBF: merged bin in TB 0, TB 1 empty, UB 2 in TBs 2-63.
    // IBF 0 (below TB 0 of the Root-IBF): UB 0 in TB 0, TB 1 empty, UB 1 in TBs 2-63.
    layout_t layout{};
    layout.top_level_max_bin_id = 0u;
    layout.max_bins.emplace_back(std::vector<size_t>{0u}, 0u);
    layout.user_bins.emplace_back(std::vector<size_t>{0u}, 0u, 1u, 0u);
    layout.user_bins.emplace_back(std::vector<size_t>{0u}, 2u, 62u, 1u);
    layout.user_bins.emplace_back(std::vector<size_t>{}, 2u, 62u, 2u);

    auto const [valid, diagnostics] = validate(layout, seqan::hibf::config{.number_of_user_bins = 3u, .tmax = 64u});
    EXPECT_TRUE(valid);
    ASSERT_EQ(diagnostics.size(), 2u);
    EXPECT_EQ(diagnostics[0].ibf, std::vector<size_t>{});
    EXPECT_EQ(diagnostics[0].message, "The Root-IBF uses technical bins 0-63, but technical bin 1 is empty.");
    EXPECT_EQ(diagnostics[1].ibf, std::vector<size_t>{0u});
    EXPECT_EQ(diagnostics[1].message, "IBF 0 uses technical bins 0-63, but technical bin 1 is empty.");
}

TEST(layout_test, validate_computed_layout)
{
    // splitmix64 finaliser
    auto scramble = [](uint64_t x)
    {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    };

    for (double const empty_bin_fraction : {0.0, 0.1})
    {
        SCOPED_TRACE(empty_bin_fraction);
        seqan::hibf::config config{.input_fn =
                                       [&](size_t const ub, seqan::hibf::insert_iterator it)
                                   {
                                       // Zipf-like sizes: The small user bins are merged, and the merged bins
                                       // contain enough user bins to be laid out hierarchically again.
                                       uint64_t const offset = static_cast<uint64_t>(ub) << 32;
                                       uint64_t const size = 20'000u / (ub + 1u) + 1u;
                                       for (uint64_t i = 0; i < size; ++i)
                                           it = scramble(offset + i);
                                   },
                                   .number_of_user_bins = 1000u,
                                   .tmax = 64u,
                                   .empty_bin_fraction = empty_bin_fraction,
                                   .disable_estimate_union = true}; // also disables rearrangement; for speed
        config.validate_and_set_defaults();

        std::vector<seqan::hibf::sketch::hyperloglog> sketches;
        std::vector<size_t> kmer_counts;
        seqan::hibf::sketch::compute_sketches(config, sketches);
        seqan::hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);

        auto const layout = seqan::hibf::layout::compute_layout(config, kmer_counts, sketches);

        EXPECT_EQ(layout.number_of_levels(), 3u);

        auto const [valid, diagnostics] = validate(layout, config);
        EXPECT_TRUE(valid);
        EXPECT_TRUE(diagnostics.empty()) << diagnostics.front(); // only evaluated on failure
    }
}
