// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides EXPECT_SAME_TYPE.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h> // for EXPECT_PRED_FORMAT2, AssertionSuccess, AssertionResult, CmpHelperEQFailure

#include <cassert>     // for assert
#include <cstddef>     // for size_t
#include <string>      // for string, allocator
#include <tuple>       // for tuple
#include <type_traits> // for type_identity

#include <hibf/test/type_name_as_string.hpp> // IWYU pragma: keep

namespace hibf::test
{
struct expect_same_type;

// https://stackoverflow.com/a/62984543
#define EXPECT_SAME_TYPE_DEPAREN(X) EXPECT_SAME_TYPE_ESC(EXPECT_SAME_TYPE_ISH X)
#define EXPECT_SAME_TYPE_ISH(...) EXPECT_SAME_TYPE_ISH __VA_ARGS__
#define EXPECT_SAME_TYPE_ESC(...) EXPECT_SAME_TYPE_ESC_(__VA_ARGS__)
#define EXPECT_SAME_TYPE_ESC_(...) EXPECT_SAME_TYPE_VAN##__VA_ARGS__
#define EXPECT_SAME_TYPE_VANEXPECT_SAME_TYPE_ISH

#define EXPECT_SAME_TYPE(val1, val2)                                                                                   \
    EXPECT_PRED_FORMAT2(::hibf::test::expect_same_type{},                                                              \
                        (std::type_identity<EXPECT_SAME_TYPE_DEPAREN(val1)>{}),                                        \
                        (std::type_identity<EXPECT_SAME_TYPE_DEPAREN(val2)>{}));

struct expect_same_type
{
    template <typename lhs_t, typename rhs_t>
    ::testing::AssertionResult operator()(std::string lhs_expression,
                                          std::string rhs_expression,
                                          std::type_identity<lhs_t>,
                                          std::type_identity<rhs_t>)
    {
        auto remove_wrap_type_identity = [](std::string str)
        {
            // EXPECT_SAME_TYPE_DEPAREN adds a space after the prefix
            std::string prefix = "std::type_identity< ";
            std::string suffix = ">{}";

            size_t str_start = str.find(prefix) + prefix.size();
            size_t str_end = str.rfind(suffix);
            assert(str_end >= str_start);

            return str.substr(str_start, str_end - str_start);
        };

        if (std::is_same_v<lhs_t, rhs_t>)
            return ::testing::AssertionSuccess();

        return ::testing::internal::CmpHelperEQFailure(remove_wrap_type_identity(lhs_expression).c_str(),
                                                       remove_wrap_type_identity(rhs_expression).c_str(),
                                                       hibf::detail::type_name_as_string<lhs_t>,
                                                       hibf::detail::type_name_as_string<rhs_t>);
    }
};

} // namespace hibf::test
