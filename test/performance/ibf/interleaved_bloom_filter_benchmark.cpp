// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for State, Benchmark, AddCustomContext, Counter, BENCHMARK

#include <algorithm>  // for __fn, generate
#include <cmath>      // for log, ceil, exp
#include <cstddef>    // for size_t
#include <functional> // for equal_to
#include <random>     // for uniform_int_distribution, mt19937_64
#include <ranges>     // for transform_view, iota_view, __range_adaptor_closure_t, __fn
#include <string>     // for to_string, basic_string
#include <tuple>      // for tuple, make_tuple
#include <utility>    // for move, pair
#include <vector>     // for vector

#include <hibf/contrib/robin_hood.hpp>              // for hash, unordered_map
#include <hibf/contrib/std/chunk_view.hpp>          // for chunk, chunk_fn, chunk_view
#include <hibf/contrib/std/detail/adaptor_base.hpp> // for operator|
#include <hibf/contrib/std/to.hpp>                  // for to
#include <hibf/interleaved_bloom_filter.hpp>        // for bin_index, interleaved_bloom_filter, bin_count, bin_size
#include <hibf/misc/divide_and_ceil.hpp>            // for divide_and_ceil
#include <hibf/platform.hpp>                        // for HIBF_HAS_AVX512
#include <hibf/test/bytes.hpp>                      // for operator""_MiB

using namespace seqan::hibf::test::literals;
static constexpr size_t total_ibf_size_in_bytes{1_MiB};
static constexpr size_t number_of_hash_functions{2u};
static constexpr double false_positive_rate{0.05};
// This computes how many elements need to be inserted into the IBF to achieve the desired false positive rate for the
// given size.
// The `number_of_elements` many generated values are used for both constructing and querying the IBF.
static /* cmath not constexpr in libc++ */ size_t number_of_elements = []()
{
    size_t const bits = 8u * total_ibf_size_in_bytes;
    double const numerator = -std::log(1 - std::exp(std::log(false_positive_rate) / number_of_hash_functions)) * bits;
    return std::ceil(numerator / number_of_hash_functions);
}();

// We cache the generated values and constructed IBFs such that they are reused for all benchmarks and multiple runs of
// the same benchmark. Google benchmark may rerun the same benchmark multiple times in order to reach the requested
// minimum benchmark time or when multiple runs are requested via the --benchmark_repetitions flag.
// Caching also means that we access the values via `const &`: `auto const & [values, ibf] = set_up(state);`.
// If a benchmark modifies the IBFs, it needs to make a copy of the IBF (clear_benchmark) or
// construct a new one (emplace_benchmark).
static robin_hood::unordered_map<size_t, std::tuple<std::vector<size_t>, seqan::hibf::interleaved_bloom_filter>>
    cache{};

auto set_up(::benchmark::State const & state)
{
    size_t const bins = state.range(0);
    size_t const bits = 8u * total_ibf_size_in_bytes / bins;

    if (auto it = cache.find(bins); it != cache.end())
        return it->second;

    std::vector<size_t> const values = []()
    {
        std::mt19937_64 engine{0ULL};
        std::uniform_int_distribution<size_t> distribution{};

        auto gen = [&]()
        {
            return distribution(engine);
        };

        std::vector<size_t> result(number_of_elements);
        std::ranges::generate(result, gen);

        return result;
    }();

    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{bins},
                                              seqan::hibf::bin_size{bits},
                                              seqan::hibf::hash_function_count{number_of_hash_functions}};

    size_t const chunk_size = seqan::hibf::divide_and_ceil(number_of_elements, bins);
    size_t bin_index = 0u;
    for (auto && chunk : seqan::stl::views::chunk(values, chunk_size))
    {
        for (auto value : chunk)
            ibf.emplace(value, seqan::hibf::bin_index{bin_index});
        ++bin_index;
    }

    auto it = cache.emplace(bins, std::make_tuple(std::move(values), std::move(ibf)));
    return it.first->second;
}

inline benchmark::Counter elements_per_second(size_t const count)
{
    return benchmark::Counter(count, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

void emplace_benchmark(::benchmark::State & state)
{
    auto const & [values, original_ibf] = set_up(state);

    size_t const bins = state.range(0);
    size_t const chunk_size = seqan::hibf::divide_and_ceil(number_of_elements, bins);

    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{original_ibf.bin_count()},
                                              seqan::hibf::bin_size{original_ibf.bin_size()},
                                              seqan::hibf::hash_function_count{original_ibf.hash_function_count()}};

    for (auto _ : state)
    {
        size_t bin_index = 0u;
        for (auto && chunk : seqan::stl::views::chunk(values, chunk_size))
        {
            for (auto value : chunk)
                ibf.emplace(value, seqan::hibf::bin_index{bin_index});
            ++bin_index;
        }
    }

    state.counters["elements"] = elements_per_second(number_of_elements);
}

void clear_benchmark(::benchmark::State & state)
{
    auto const & [values, original_ibf] = set_up(state);
    (void)values;

    size_t const bins = state.range(0);

    seqan::hibf::interleaved_bloom_filter ibf{original_ibf};

    for (auto _ : state)
    {
        for (size_t i = 0; i < bins; ++i)
            ibf.clear(seqan::hibf::bin_index{i});
    }

    state.counters["bins"] = elements_per_second(bins);
}

void clear_range_benchmark(::benchmark::State & state)
{
    auto const & [values, original_ibf] = set_up(state);
    (void)values;

    size_t const bins = state.range(0);

    seqan::hibf::interleaved_bloom_filter ibf{original_ibf};

    std::vector<seqan::hibf::bin_index> bin_range = std::views::iota(0u, static_cast<size_t>(state.range(0)))
                                                  | std::views::transform(
                                                        [](size_t i)
                                                        {
                                                            return seqan::hibf::bin_index{i};
                                                        })
                                                  | seqan::stl::ranges::to<std::vector>();

    for (auto _ : state)
    {
        ibf.clear(bin_range);
    }

    state.counters["bins"] = elements_per_second(bins);
}

void bulk_contains_benchmark(::benchmark::State & state)
{
    auto const & [values, ibf] = set_up(state);

    auto agent = ibf.membership_agent();
    for (auto _ : state)
    {
        for (auto hash : values)
        {
            [[maybe_unused]] auto & res = agent.bulk_contains(hash);
            benchmark::ClobberMemory();
        }
    }

    state.counters["elements"] = elements_per_second(number_of_elements);
}

void bulk_count_benchmark(::benchmark::State & state)
{
    auto const & [values, ibf] = set_up(state);

    auto agent = ibf.counting_agent();
    for (auto _ : state)
    {
        [[maybe_unused]] auto & res = agent.bulk_count(values);
        benchmark::ClobberMemory();
    }

    state.counters["elements"] = elements_per_second(number_of_elements);
}

BENCHMARK(emplace_benchmark)->RangeMultiplier(2)->Range(64, 1024);
BENCHMARK(clear_benchmark)->RangeMultiplier(2)->Range(64, 1024);
BENCHMARK(clear_range_benchmark)->RangeMultiplier(2)->Range(64, 1024);
BENCHMARK(bulk_contains_benchmark)->RangeMultiplier(2)->Range(64, 1024);
BENCHMARK(bulk_count_benchmark)->RangeMultiplier(2)->Range(64, 1024);

// This is a hack to add custom context information to the benchmark output.
// The alternative would be to do it in the main(). However, this would require
// not using the BENCHMARK_MAIN macro.
[[maybe_unused]] static bool foo = []()
{
    benchmark::AddCustomContext("IBF size in bytes", std::to_string(total_ibf_size_in_bytes));
    benchmark::AddCustomContext("Number of hash functions", std::to_string(number_of_hash_functions));
    benchmark::AddCustomContext("False positive rate", std::to_string(false_positive_rate));
    benchmark::AddCustomContext("Number of elements", std::to_string(number_of_elements));
    benchmark::AddCustomContext("HIBF_HAS_AVX512", HIBF_HAS_AVX512 ? "true" : "false");
    benchmark::AddCustomContext("AVX512 support",
#if __AVX512F__ && __AVX512BW__
                                "true");
#else
                                "false");
#endif
    return true;
}();

BENCHMARK_MAIN();
