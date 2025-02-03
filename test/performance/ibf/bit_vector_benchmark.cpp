// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

// Authored by: Rene Rahn <rene.rahn AT fu-berlin.de>

#include <benchmark/benchmark.h> // for Benchmark, BENCHMARK_CAPTURE, DoNotOptimize, State, BENCHMARK_MAIN

#include <algorithm>   // for __fn, generate
#include <cstddef>     // for size_t
#include <cstdint>     // for int32_t, uint8_t
#include <random>      // for uniform_int_distribution, mt19937_64
#include <string>      // for basic_string
#include <type_traits> // for invoke_result_t
#include <utility>     // for move, pair

#include <hibf/misc/bit_vector.hpp> // for bit_vector

seqan::hibf::bit_vector generate_bit_vector(size_t const size, size_t const seed)
{
    std::mt19937_64 engine{seed};
    std::uniform_int_distribution<uint8_t> dist{0u, 1u};

    auto gen = [&dist, &engine]()
    {
        return dist(engine);
    };
    seqan::hibf::bit_vector vec(size);
    std::ranges::generate(vec, gen);

    return vec;
}

auto generate_bit_vector_pair(size_t const size)
{
    seqan::hibf::bit_vector random_bit_vector_first = generate_bit_vector(size, 0u);
    seqan::hibf::bit_vector random_bit_vector_second = generate_bit_vector(size, size);

    return std::pair{std::move(random_bit_vector_first), std::move(random_bit_vector_second)};
}

template <typename operation_t>
void random_bit_vector(benchmark::State & state, operation_t operation)
{
    auto [test_vector_lhs, test_vector_rhs] = generate_bit_vector_pair(state.range(0));

    using result_t = std::invoke_result_t<operation_t, seqan::hibf::bit_vector &, seqan::hibf::bit_vector &>;

    for (auto _ : state)
    {
        result_t result = operation(test_vector_lhs, test_vector_rhs);
        benchmark::DoNotOptimize(result);
    }
}

template <typename operation_t>
void one_set_bit_vector(benchmark::State & state, operation_t operation)
{
    seqan::hibf::bit_vector vec(state.range(0));
    vec.back() = true;

    for (auto _ : state)
    {
        bool result = operation(vec);
        benchmark::DoNotOptimize(result);
    }
}

// ----------------------------------------------------------------------------
// Benchmark operations
// ----------------------------------------------------------------------------

inline auto binary_and = []<typename bv_t>(bv_t & lhs, bv_t & rhs) constexpr -> bv_t &
{
    return lhs &= rhs;
};

inline auto binary_or = []<typename bv_t>(bv_t & lhs, bv_t & rhs) constexpr -> bv_t &
{
    return lhs |= rhs;
};

inline auto binary_xor = []<typename bv_t>(bv_t & lhs, bv_t & rhs) constexpr -> bv_t &
{
    return lhs ^= rhs;
};

inline auto binary_not = []<typename bv_t>(bv_t & lhs, [[maybe_unused]] auto const &... args) constexpr -> bv_t
{
    return ~lhs;
};

inline auto binary_flip = []<typename bv_t>(bv_t & lhs, [[maybe_unused]] auto const &... args) constexpr -> bv_t &
{
    return lhs.flip();
};

inline auto none_fn = [](auto & lhs) constexpr -> bool
{
    return lhs.none();
};

inline auto all_fn = [](auto & lhs) constexpr -> bool
{
    return lhs.all();
};

inline auto any_fn = [](auto & lhs) constexpr -> bool
{
    return lhs.any();
};

#if 0
static constexpr int32_t min_range = 262'144;
static constexpr int32_t max_range = 1'048'576;
#else
static constexpr int32_t min_range = 1'024;
static constexpr int32_t max_range = 2'048;
#endif
static constexpr int32_t range_multiplier = 2;

// ----------------------------------------------------------------------------
// Benchmark
// ----------------------------------------------------------------------------

BENCHMARK_CAPTURE(random_bit_vector, and, binary_and)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);
BENCHMARK_CAPTURE(random_bit_vector, or, binary_or)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);
BENCHMARK_CAPTURE(random_bit_vector, xor, binary_xor)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);
BENCHMARK_CAPTURE(random_bit_vector, not, binary_not)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);
BENCHMARK_CAPTURE(random_bit_vector, flip, binary_flip)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

BENCHMARK_CAPTURE(one_set_bit_vector, none, none_fn)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);
BENCHMARK_CAPTURE(one_set_bit_vector, all, all_fn)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);
BENCHMARK_CAPTURE(one_set_bit_vector, any, any_fn)->RangeMultiplier(range_multiplier)->Range(min_range, max_range);

// ----------------------------------------------------------------------------
// Run benchmark
// ----------------------------------------------------------------------------

BENCHMARK_MAIN();
