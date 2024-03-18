// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for State, DoNotOptimize, Benchmark, BENCHMARK, BENCHMARK_MAIN

#include <cinttypes> // for uint8_t

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
BENCHMARK(merge)->DenseRange(5, 12, 1);
BENCHMARK(estimate)->DenseRange(5, 12, 1);
BENCHMARK(merge_and_estimate)->DenseRange(5, 12, 1);

BENCHMARK_MAIN();
