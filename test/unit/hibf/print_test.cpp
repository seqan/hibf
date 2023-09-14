#include <gtest/gtest.h>

#include <hibf/misc/print.hpp>

template <typename t>
using print_test = ::testing::Test;

using test_types = ::testing::Types<seqan::hibf::binning_bitvector,
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
    if constexpr (std::same_as<TypeParam, seqan::hibf::binning_bitvector>)
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
    if constexpr (std::same_as<TypeParam, seqan::hibf::binning_bitvector>)
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
