// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides test utilities for std::ranges::range types.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h> // for EXPECT_PRED_FORMAT2, AssertionSuccess, AssertionResult, CmpHelperEQFailure

#include <algorithm> // for copy, equal
#include <concepts>  // for same_as
#include <iterator>  // for back_inserter
#include <ostream>   // for operator<<, ostream
#include <ranges>    // for range, range_value_t
#include <utility>   // for forward
#include <vector>    // for vector

#include <hibf/platform.hpp>

namespace seqan::hibf::test
{
struct expect_range_eq;
} // namespace seqan::hibf::test

namespace std
{

template <typename t>
    requires std::same_as<t, std::vector<char>>
void PrintTo(t const & value, std::ostream * out)
{
    for (char c : value)
        *out << c;
}

} // namespace std

namespace seqan::hibf::test
{

#define EXPECT_RANGE_EQ(val1, val2) EXPECT_PRED_FORMAT2(::seqan::hibf::test::expect_range_eq{}, val1, val2);

struct expect_range_eq
{
    template <typename rng_t>
    auto copy_range(rng_t && rng)
    {
        using value_t = std::ranges::range_value_t<rng_t>;
        std::vector<value_t> rng_copy{};
        std::ranges::copy(rng, std::back_inserter(rng_copy));
        return rng_copy;
    }

    template <std::ranges::range lhs_t, std::ranges::range rhs_t>
    ::testing::AssertionResult
    operator()(char const * lhs_expression, char const * rhs_expression, lhs_t && lhs, rhs_t && rhs)
    {
        std::vector lhs_copy = copy_range(std::forward<lhs_t>(lhs));
        std::vector rhs_copy = copy_range(std::forward<rhs_t>(rhs));

        if (std::ranges::equal(lhs_copy, rhs_copy))
            return ::testing::AssertionSuccess();

        return ::testing::internal::CmpHelperEQFailure(lhs_expression, rhs_expression, lhs_copy, rhs_copy);
    }
};

} // namespace seqan::hibf::test
