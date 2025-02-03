// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HIBF_HAS_AVX512
#define HIBF_HAS_AVX512 0
#endif

#include <gtest/gtest.h> // for Message, TYPED_TEST, TestPartResult, ASSERT_EQ, EXPECT_EQ, Types

#include <algorithm> // for fill, fill_n
#include <concepts>  // for signed_integral
#include <cstddef>   // for size_t
#include <cstdint>   // for int16_t, int32_t, int64_t, int8_t, uint16_t, uint32_t, uint64_t
#include <limits>    // for numeric_limits
#include <string>    // for basic_string

#include <hibf/misc/bit_vector.hpp>      // for bit_vector
#include <hibf/misc/counting_vector.hpp> // for counting_vector
#include <hibf/platform.hpp>             // for HIBF_COMPILER_IS_GCC, _LIBCPP_HAS_NO_ASAN, _LIBCPP_VERSION

template <typename TypeParam>
class counting_vector_test : public ::testing::Test
{
protected:
    seqan::hibf::counting_vector<TypeParam> counting_vector;
    seqan::hibf::bit_vector bit_vector;
    seqan::hibf::counting_vector<TypeParam> expected;

    void SetUp() override
    {
        counting_vector.resize(2048);
        bit_vector.resize(2048);
        expected.resize(2048);
    }

    void check_add()
    {
        counting_vector += bit_vector;
        ASSERT_EQ(counting_vector.size(), expected.size());
        EXPECT_EQ(counting_vector, expected);
    }

    void check_sub()
    {
        counting_vector -= bit_vector;
        ASSERT_EQ(counting_vector.size(), expected.size());
        EXPECT_EQ(counting_vector, expected);
    }

    void annotate_llvm_asan() const
    {
#if defined(_LIBCPP_VERSION) && !defined(_LIBCPP_HAS_NO_ASAN)
        __sanitizer_annotate_contiguous_container(counting_vector.data(),
                                                  counting_vector.data() + counting_vector.capacity(),
                                                  counting_vector.data() + counting_vector.size(),
                                                  counting_vector.data() + counting_vector.capacity());
#endif
    }
};

using test_types = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t, int8_t, int16_t, int32_t, int64_t>;

TYPED_TEST_SUITE(counting_vector_test, test_types);

TYPED_TEST(counting_vector_test, empty)
{
    this->counting_vector.clear();
    this->bit_vector.clear();
    this->expected.clear();

    this->check_add();
    this->check_sub();
}

TYPED_TEST(counting_vector_test, no_bits_set)
{
    this->check_add();
    this->check_sub();
}

TYPED_TEST(counting_vector_test, all_bits_set)
{
    constexpr TypeParam max_value = std::numeric_limits<TypeParam>::max();

    std::ranges::fill(this->counting_vector, max_value - 1u);
    std::ranges::fill(this->bit_vector, true);

    std::ranges::fill(this->expected, max_value);
    this->check_add();

    std::ranges::fill(this->expected, max_value - 1u);
    this->check_sub();

    if constexpr (std::signed_integral<TypeParam>)
    {
        std::ranges::fill(this->counting_vector, 0);

        std::ranges::fill(this->expected, -1);
        this->check_sub();
    }
}

TYPED_TEST(counting_vector_test, every_other_bit_set)
{
    for (size_t i = 0; i < this->bit_vector.size(); ++i)
        this->bit_vector[i] = i % 2;

    for (size_t i = 0; i < this->expected.size(); ++i)
        this->expected[i] = i % 2;
    this->check_add();

    std::ranges::fill(this->expected, TypeParam{});
    this->check_sub();

    if constexpr (std::signed_integral<TypeParam>)
    {
        std::ranges::fill(this->counting_vector, 0);

        for (size_t i = 0; i < this->expected.size(); ++i)
            this->expected[i] = -(i % 2);
        this->check_sub();
    }
}

TYPED_TEST(counting_vector_test, first_and_last_set)
{
    this->bit_vector.front() = this->bit_vector.back() = true;

    ++this->expected.front();
    ++this->expected.back();
    this->check_add();

    --this->expected.front();
    --this->expected.back();
    this->check_sub();

    if constexpr (std::signed_integral<TypeParam>)
    {
        --this->expected.front();
        --this->expected.back();
        this->check_sub();
    }
}

TYPED_TEST(counting_vector_test, counting_vector_is_bigger)
{
    constexpr TypeParam max_value = std::numeric_limits<TypeParam>::max();

    this->counting_vector.resize(4096, max_value);
    std::ranges::fill(this->bit_vector, true);

    this->expected.resize(4096, max_value);
    std::ranges::fill_n(this->expected.begin(), this->expected.size() / 2u, TypeParam{1u});
    this->check_add();

    std::ranges::fill_n(this->expected.begin(), this->expected.size() / 2u, TypeParam{});
    this->check_sub();
}
TYPED_TEST(counting_vector_test, size_not_divisible_by_64)
{
    constexpr TypeParam max_value = std::numeric_limits<TypeParam>::max();

    this->counting_vector.resize(2033);
    std::ranges::fill(this->counting_vector, max_value - 1u);
    this->annotate_llvm_asan();

    this->bit_vector.resize(2033);
    std::ranges::fill(this->bit_vector, true);

    this->expected.resize(2033);
    std::ranges::fill(this->expected, max_value);
    this->check_add();

    std::ranges::fill(this->expected, max_value - 1u);
    this->check_sub();
}

#if !(HIBF_COMPILER_IS_GCC && defined(__SANITIZE_ADDRESS__))
TYPED_TEST(counting_vector_test, overflow)
{
    constexpr TypeParam max_value = std::numeric_limits<TypeParam>::max();

    std::ranges::fill(this->counting_vector, max_value);
    this->counting_vector.resize(2033);
    this->annotate_llvm_asan();

    std::ranges::fill(this->bit_vector, true);
    this->bit_vector.resize(2033);
    std::ranges::fill(this->bit_vector, false);

    this->expected.resize(2033);
    std::ranges::fill(this->expected, max_value);
    this->check_add(); // counting_vector[2033...2047] may overflow but it does not affect the actual results

    std::ranges::fill(this->expected, max_value);
    this->check_sub(); // counting_vector[2033...2047] may underflow but it does not affect the actual results
}
#endif
