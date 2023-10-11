#include <gtest/gtest.h> // for Test, Message, TestPartResult, EXPECT_EQ, TestInfo

#include <cinttypes>   // for uint8_t
#include <cstddef>     // for size_t
#include <random>      // for uniform_int_distribution, mt19937_64
#include <string>      // for basic_string, string
#include <string_view> // for string_view
#include <vector>      // for allocator, vector

#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog

TEST(estimate_kmer_counts_test, small_example)
{
    seqan::hibf::sketch::hyperloglog sketch(12);

    std::uniform_int_distribution<uint64_t> distribution{};
    std::mt19937_64 engine{0u};

    for (size_t i = 0; i < 1500; ++i)
        sketch.add(distribution(engine));

    std::vector<seqan::hibf::sketch::hyperloglog> sketches{sketch, sketch};
    std::vector<size_t> kmer_counts;

    seqan::hibf::sketch::estimate_kmer_counts(sketches, kmer_counts);

    ASSERT_EQ(kmer_counts.size(), 2);
    EXPECT_EQ(kmer_counts[0], 1498);
    EXPECT_EQ(kmer_counts[1], 1498);
}
