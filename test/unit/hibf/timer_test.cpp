// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Message, AssertionResult, TestPartResult, CmpHelperGE, EXPECT_TRUE, CmpHelper...

#include <chrono>      // for operator""ms
#include <string>      // for basic_string
#include <thread>      // for sleep_for
#include <type_traits> // for is_copy_assignable_v, is_copy_constructible_v, is_default_constructible_v

#include <hibf/misc/timer.hpp> // for concurrent_timer, serial_timer

static inline void waste_time()
{
    using namespace std::chrono_literals;
    std::this_thread::sleep_for(1ms);
}

template <typename t>
using timer_test = ::testing::Test;

using test_types = ::testing::Types<seqan::hibf::serial_timer, seqan::hibf::concurrent_timer>;

TYPED_TEST_SUITE(timer_test, test_types);

TYPED_TEST(timer_test, concepts)
{
    EXPECT_TRUE(std::is_default_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_copy_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_move_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_copy_assignable_v<TypeParam>);
    EXPECT_TRUE(std::is_move_assignable_v<TypeParam>);
    EXPECT_TRUE(std::is_destructible_v<TypeParam>);
}

#ifndef NDEBUG
TYPED_TEST(timer_test, assertions)
{
    TypeParam timer{};
    EXPECT_DEATH(timer.stop(), "stop_point >= start_point");
    EXPECT_DEATH(timer.avg_in_seconds(), "count(\\.load\\(\\))? > 0u");
}
#endif

TYPED_TEST(timer_test, in_seconds)
{
    TypeParam timer{};
    EXPECT_EQ(timer.in_seconds(), 0.0);
    timer.start();
    waste_time();
    timer.stop();
    double const first_time = timer.in_seconds();
    EXPECT_GE(first_time, 0.001);

    timer.start();
    waste_time();
    timer.stop();
    double const second_time = timer.in_seconds();
    EXPECT_GE(second_time, 0.002);
    EXPECT_GT(second_time, first_time);
}

TYPED_TEST(timer_test, max_in_seconds)
{
    TypeParam timer{};
    timer.start();
    timer.stop();
    double const first_max = timer.max_in_seconds();
    EXPECT_GE(first_max, 0.0);
    EXPECT_EQ(first_max, timer.in_seconds());

    timer.start();
    waste_time();
    timer.stop();
    double const second_max = timer.max_in_seconds();
    EXPECT_GT(second_max, first_max);
}

TYPED_TEST(timer_test, avg_in_seconds)
{
    TypeParam timer{};
    timer.start();
    timer.stop();
    double const first_avg = timer.avg_in_seconds();
    EXPECT_GE(first_avg, 0.0);
    EXPECT_EQ(first_avg, timer.in_seconds());

    timer.start();
    waste_time();
    timer.stop();
    double const second_avg = timer.avg_in_seconds();
    EXPECT_EQ(second_avg, timer.in_seconds() / 2.0);
    EXPECT_GT(second_avg, first_avg);
}

TYPED_TEST(timer_test, add)
{
    TypeParam timer{};

    seqan::hibf::serial_timer timer2{};
    timer2.start();
    waste_time();
    timer2.stop();

    seqan::hibf::concurrent_timer timer3{};
    timer3.start();
    waste_time();
    timer3.stop();

    EXPECT_EQ(timer.in_seconds(), 0.0);

    timer += timer2;
    EXPECT_GE(timer.in_seconds(), 0.001);

    timer += timer3;
    EXPECT_GE(timer.in_seconds(), 0.002);
}

TYPED_TEST(timer_test, comparison)
{
    TypeParam timer1{};
    seqan::hibf::serial_timer timer2{};
    seqan::hibf::concurrent_timer timer3{};
    EXPECT_TRUE(timer1 == timer2);
    EXPECT_TRUE(timer1 == timer3);
}
