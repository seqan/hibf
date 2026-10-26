// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for Benchmark, DoNotOptimize, ClobberMemory, State, BENCHMARK, BENCHMARK_MAIN

#include <cstddef> // for size_t
#include <cstdint> // for uint64_t, uint8_t
#include <random>  // for mt19937_64, uniform_int_distribution
#include <vector>  // for vector

#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

static void add(benchmark::State & state)
{
    uint8_t const bits = static_cast<uint8_t>(state.range(0));
    seqan::hibf::sketch::hyperloglog hll{bits};
    for (auto _ : state)
    {
        hll.add(23);
    }
    benchmark::DoNotOptimize(hll);
}

static void add_many(benchmark::State & state)
{
    static constexpr size_t number_of_values = 1ULL << 20;
    static std::vector<uint64_t> const values = []()
    {
        std::mt19937_64 engine{42u};
        std::uniform_int_distribution<uint64_t> distribution{};
        std::vector<uint64_t> result(number_of_values);
        for (uint64_t & value : result)
            value = distribution(engine);
        return result;
    }();

    uint8_t const bits = static_cast<uint8_t>(state.range(0));
    seqan::hibf::sketch::hyperloglog hll{bits};
    for (auto _ : state)
    {
        for (uint64_t const value : values)
            hll.add(value);
        benchmark::DoNotOptimize(hll);
        benchmark::ClobberMemory();
    }

    state.SetItemsProcessed(state.iterations() * number_of_values);
}

static void merge(benchmark::State & state)
{
    uint8_t const bits = static_cast<uint8_t>(state.range(0));
    seqan::hibf::sketch::hyperloglog hll{bits};
    for (auto _ : state)
    {
        hll.merge(hll);
    }
    benchmark::DoNotOptimize(hll);
}

static void estimate(benchmark::State & state)
{
    uint8_t const bits = static_cast<uint8_t>(state.range(0));
    seqan::hibf::sketch::hyperloglog hll{bits};
    double estimate{};
    for (auto _ : state)
    {
        estimate = hll.estimate();
    }
    benchmark::DoNotOptimize(hll);
    benchmark::DoNotOptimize(estimate);
}

static void merge_and_estimate(benchmark::State & state)
{
    uint8_t const bits = static_cast<uint8_t>(state.range(0));
    seqan::hibf::sketch::hyperloglog hll{bits};
    double estimate{};
    for (auto _ : state)
    {
        estimate = hll.merge_and_estimate(hll);
    }
    benchmark::DoNotOptimize(hll);
    benchmark::DoNotOptimize(estimate);
}

BENCHMARK(add)->DenseRange(5, 12, 1);
BENCHMARK(add_many)->DenseRange(5, 12, 1);
BENCHMARK(merge)->DenseRange(5, 12, 1);
BENCHMARK(estimate)->DenseRange(5, 12, 1);
BENCHMARK(merge_and_estimate)->DenseRange(5, 12, 1);

BENCHMARK_MAIN();
