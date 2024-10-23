// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for BENCHMARK_TEMPLATE, Benchmark, State

#include <cinttypes> // for uint8_t
#include <cstddef>   // for size_t
#include <vector>    // for vector

#include <hibf/config.hpp>                  // for config, insert_iterator
#include <hibf/sketch/compute_sketches.hpp> // for compute_sketches
#include <hibf/sketch/hyperloglog.hpp>      // for hyperloglog
#include <hibf/sketch/minhashes.hpp>        // for minhashes

inline benchmark::Counter elements_per_second(size_t const count)
{
    return benchmark::Counter(count, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

enum class sketch : uint8_t
{
    Hyperloglog,
    MinHashes
};

template <sketch sketch_t>
void compute_sketches(benchmark::State & state)
{
    static constexpr uint64_t elements_per_bin = 10000;
    auto create_hashes = [&](size_t const ub_id, seqan::hibf::insert_iterator it)
    {
        // 0 = [0, 10000]
        // 1 = [10000, 20000]
        // 1 = [20000, 30000]
        for (size_t i = ub_id * elements_per_bin; i < (ub_id + 1) * elements_per_bin; ++i)
            it = i;
    };

    seqan::hibf::config config{};
    config.number_of_user_bins = 64;
    config.input_fn = create_hashes;
    config.sketch_bits = 12;

    [[maybe_unused]] std::vector<seqan::hibf::sketch::minhashes> minhash_sketches;
    std::vector<seqan::hibf::sketch::hyperloglog> hyperloglog_sketches(config.number_of_user_bins, config.sketch_bits);

    for (auto _ : state)
    {
        if constexpr (sketch_t == sketch::MinHashes)
            seqan::hibf::sketch::compute_sketches(config, hyperloglog_sketches, minhash_sketches);
        else
            seqan::hibf::sketch::compute_sketches(config, hyperloglog_sketches);

        benchmark::DoNotOptimize(hyperloglog_sketches);
        benchmark::ClobberMemory();
    }

    state.counters["elements"] = elements_per_second(elements_per_bin * config.number_of_user_bins);
}

BENCHMARK_TEMPLATE(compute_sketches, sketch::Hyperloglog);
BENCHMARK_TEMPLATE(compute_sketches, sketch::MinHashes);
