// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for State, Benchmark, Counter, ClobberMemory, BENCHMARK

#include <algorithm>  // for copy, fill_n, __generate_fn, generate
#include <cstddef>    // for size_t
#include <functional> // for function
#include <random>     // for uniform_int_distribution, mt19937_64
#include <ranges>     // for _Size, size
#include <span>       // for span
#include <stdexcept>  // for invalid_argument
#include <string>     // for basic_string
#include <tuple>      // for make_tuple
#include <vector>     // for allocator, vector

#include <hibf/config.hpp>                                // for insert_iterator, config
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

inline benchmark::Counter hashes_per_second(size_t const count)
{
    return benchmark::Counter(count, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

static constexpr size_t hash_num{2u};
static constexpr double fpr{0.05};

#ifdef NDEBUG                                        // Release
static constexpr size_t sequence_length{1ULL << 18}; // Distributed evenly across all user bins.
#else                                                // Debug
static constexpr size_t sequence_length{1ULL << 12}; // Distributed evenly across all user bins.
#endif

auto set_up(::benchmark::State const & state)
{
    size_t const num_ub = state.range(0);
    if (sequence_length < num_ub)
        throw std::invalid_argument{"sequence_length must be >= num_ub."};
    if (sequence_length % num_ub != 0u)
        throw std::invalid_argument{"sequence_length % num_ub must be 0."};
    size_t const values_per_ub = sequence_length / num_ub;

    // Generate random values for insertion and query.
    std::vector<size_t> values(sequence_length);
    auto generator = []()
    {
        std::uniform_int_distribution<size_t> distr{0u};
        std::mt19937_64 engine{0ULL};
        return distr(engine);
    };
    std::ranges::generate(values, generator);

    auto distribute_hashes_across_ub = [&](size_t const ub_id, seqan::hibf::insert_iterator it)
    {
        for (auto value : std::span(values.begin() + ub_id * values_per_ub, values_per_ub))
            it = value;
    };

    seqan::hibf::config config{.input_fn = distribute_hashes_across_ub,
                               .number_of_user_bins = num_ub,
                               .number_of_hash_functions = hash_num,
                               .maximum_fpr = fpr,
                               .threads = 4u, // Only applies to layout and build
                               .disable_estimate_union = true};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

    return std::make_tuple(values, hibf);
}

void membership_for_benchmark_single_query(::benchmark::State & state)
{
    auto && [values, hibf] = set_up(state);
    auto agent = hibf.membership_agent();

    for (auto _ : state)
    {
        [[maybe_unused]] auto & res = agent.membership_for(values, 1u);
        benchmark::ClobberMemory();
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(values));
}

void membership_for_benchmark_multiple_queries(::benchmark::State & state)
{
    auto && [values, hibf] = set_up(state);
    auto agent = hibf.membership_agent();

    size_t const num_queries = 64u;
    size_t const values_per_query = sequence_length / num_queries;

    for (auto _ : state)
    {
        for (size_t i = 0u; i < num_queries; ++i)
        {
            auto query = std::span(values.begin() + i * values_per_query, values_per_query);
            [[maybe_unused]] auto & res = agent.membership_for(query, 1u);
            benchmark::ClobberMemory();
        }
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(values));
}

BENCHMARK(membership_for_benchmark_single_query)->RangeMultiplier(2)->Range(64, 256)->Iterations(1);
BENCHMARK(membership_for_benchmark_multiple_queries)->RangeMultiplier(2)->Range(64, 256)->Iterations(1);

BENCHMARK_MAIN();
