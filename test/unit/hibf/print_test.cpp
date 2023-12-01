// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for EXPECT_EQ, GetCapturedStderr, GetCapturedStdout, Message, TestPartR...

#include <cinttypes> // for int64_t, int16_t, int32_t, int8_t, uint16_t, uint32_t, uint64_t
#include <concepts>  // for same_as, unsigned_integral
#include <iostream>  // for cerr
#include <ranges>    // for range_value_t
#include <vector>    // for allocator, vector

#include <hibf/misc/bit_vector.hpp>      // for bit_vector
#include <hibf/misc/counting_vector.hpp> // for counting_vector
#include <hibf/misc/print.hpp>           // for print

template <typename t>
using print_test = ::testing::Test;

using test_types = ::testing::Types<seqan::hibf::bit_vector,
                                    seqan::hibf::counting_vector<uint8_t>,
                                    seqan::hibf::counting_vector<uint16_t>,
                                    seqan::hibf::counting_vector<uint32_t>,
                                    seqan::hibf::counting_vector<uint64_t>,
                                    seqan::hibf::counting_vector<int8_t>,
                                    seqan::hibf::counting_vector<int16_t>,
                                    seqan::hibf::counting_vector<int32_t>,
                                    seqan::hibf::counting_vector<int64_t>,
                                    std::vector<int64_t>>;

TYPED_TEST_SUITE(print_test, test_types);

TYPED_TEST(print_test, empty)
{
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    seqan::hibf::print(TypeParam{});
    EXPECT_EQ((testing::internal::GetCapturedStdout()), "[]\n");
    EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
}

TYPED_TEST(print_test, to_stdout)
{
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    if constexpr (std::same_as<TypeParam, seqan::hibf::bit_vector>)
    {
        TypeParam vector(5u);
        vector[0] = vector[2] = vector[4] = true;
        seqan::hibf::print(vector);
        EXPECT_EQ((testing::internal::GetCapturedStdout()), "[1,0,1,0,1]\n");
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }
    else if constexpr (std::unsigned_integral<std::ranges::range_value_t<TypeParam>>)
    {
        TypeParam vector{10u, 5u, 100u, 0u, 20u};
        seqan::hibf::print(vector);
        EXPECT_EQ((testing::internal::GetCapturedStdout()), "[10,5,100,0,20]\n");
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }
    else
    {
        TypeParam vector{-100, 2, -3, 0, 125};
        seqan::hibf::print(vector);
        EXPECT_EQ((testing::internal::GetCapturedStdout()), "[-100,2,-3,0,125]\n");
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "");
    }
}

TYPED_TEST(print_test, to_stderr)
{
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    if constexpr (std::same_as<TypeParam, seqan::hibf::bit_vector>)
    {
        TypeParam vector(5u);
        vector[0] = vector[2] = vector[4] = true;
        seqan::hibf::print(vector, std::cerr);
        EXPECT_EQ((testing::internal::GetCapturedStdout()), "");
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "[1,0,1,0,1]\n");
    }
    else if constexpr (std::unsigned_integral<std::ranges::range_value_t<TypeParam>>)
    {
        TypeParam vector{10u, 5u, 100u, 0u, 20u};
        seqan::hibf::print(vector, std::cerr);
        EXPECT_EQ((testing::internal::GetCapturedStdout()), "");
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "[10,5,100,0,20]\n");
    }
    else
    {
        TypeParam vector{-100, 2, -3, 0, 125};
        seqan::hibf::print(vector, std::cerr);
        EXPECT_EQ((testing::internal::GetCapturedStdout()), "");
        EXPECT_EQ((testing::internal::GetCapturedStderr()), "[-100,2,-3,0,125]\n");
    }
}
