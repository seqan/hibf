// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for State, Benchmark, Counter, BENCHMARK, ClobberMemory, BEN...

#include <algorithm> // for copy, fill_n, __generate_fn, generate
#include <cstddef>   // for size_t
#include <limits>    // for numeric_limits
#include <random>    // for uniform_int_distribution, mt19937_64
#include <ranges>    // for transform_view, all_t, iota_view, _Partial, _Transform
#include <string>    // for basic_string
#include <tuple>     // for make_tuple
#include <vector>    // for vector, allocator

#include <hibf/contrib/std/detail/adaptor_base.hpp> // for operator|
#include <hibf/contrib/std/pair.hpp>                // for operator==, pair
#include <hibf/contrib/std/to.hpp>                  // for to
#include <hibf/contrib/std/zip_view.hpp>            // for zip_view, operator==, zip, zip_fn
#include <hibf/interleaved_bloom_filter.hpp>        // for bin_index, interleaved_bloom_filter, bin_count, bin_size

inline benchmark::Counter hashes_per_second(size_t const count)
{
    return benchmark::Counter(count, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

#if 1
static void arguments(benchmark::internal::Benchmark * b)
{
    // Total size: 1MiB
    // bins, bin_size, hash_num, sequence_length
    b->Args({64, 1LL << 17, 2, 1LL << 17});
    b->Args({128, 1LL << 16, 2, 1LL << 17});
    b->Args({192, 1LL << 16, 2, 1LL << 17});
    b->Args({256, 1LL << 15, 2, 1LL << 17});
    // b->Args({320, 1LL << 15, 2, 1LL << 17});
    // b->Args({384, 1LL << 14, 2, 1LL << 17});
    // b->Args({448, 1LL << 14, 2, 1LL << 17});
    // b->Args({512, 1LL << 13, 2, 1LL << 17});
    b->Args({1024, 1LL << 10, 2, 1LL << 17});
}
#else
static void arguments(benchmark::internal::Benchmark * b)
{
    // Total size: 1GiB
    // bins, bin_size, hash_num, sequence_length
    b->Args({64, 1LL << 27, 2, 1LL << 27});
    b->Args({128, 1LL << 26, 2, 1LL << 27});
    b->Args({192, 1LL << 26, 2, 1LL << 27});
    b->Args({256, 1LL << 25, 2, 1LL << 27});
    b->Args({1024, 1LL << 20, 2, 1LL << 27});
}
#endif

auto set_up(::benchmark::State const & state)
{
    size_t const bins = state.range(0);
    size_t const bits = state.range(1);
    size_t const hash_num = state.range(2);
    size_t const sequence_length = state.range(3);

    auto generate = [sequence_length](size_t const max_value = std::numeric_limits<size_t>::max())
    {
        auto generator = [max_value]()
        {
            std::uniform_int_distribution<size_t> distr{0u, max_value};
            std::mt19937_64 engine{0ULL};
            return distr(engine);
        };
        std::vector<size_t> result(sequence_length);

        std::ranges::generate(result, generator);
        return result;
    };

    std::vector<size_t> const bin_indices{generate(bins - 1)};
    std::vector<size_t> const hash_values{generate()};

    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{bins},
                                              seqan::hibf::bin_size{bits},
                                              seqan::hibf::hash_function_count{hash_num}};

    return std::make_tuple(bin_indices, hash_values, ibf);
}

void emplace_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] = set_up(state);

    for (auto _ : state)
    {
        for (auto [hash, bin] : seqan::stl::views::zip(hash_values, bin_indices))
            ibf.emplace(hash, seqan::hibf::bin_index{bin});
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

void clear_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] = set_up(state);
    (void)bin_indices;
    (void)hash_values;

    std::vector<seqan::hibf::bin_index> bin_range = std::views::iota(0u, static_cast<size_t>(state.range(0)))
                                                  | std::views::transform(
                                                        [](size_t i)
                                                        {
                                                            return seqan::hibf::bin_index{i};
                                                        })
                                                  | seqan::stl::ranges::to<std::vector>();

    for (auto _ : state)
    {
        for (auto bin : bin_range)
            ibf.clear(bin);
    }

    state.counters["bins/sec"] = hashes_per_second(std::ranges::size(bin_range));
}

void clear_range_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] = set_up(state);
    (void)bin_indices;
    (void)hash_values;

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

    state.counters["bins/sec"] = hashes_per_second(std::ranges::size(bin_range));
}

void bulk_contains_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] = set_up(state);

    for (auto [hash, bin] : seqan::stl::views::zip(hash_values, bin_indices))
        ibf.emplace(hash, seqan::hibf::bin_index{bin});

    auto agent = ibf.membership_agent();
    for (auto _ : state)
    {
        for (auto hash : hash_values)
        {
            [[maybe_unused]] auto & res = agent.bulk_contains(hash);
            benchmark::ClobberMemory();
        }
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

void bulk_count_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] = set_up(state);

    for (auto [hash, bin] : seqan::stl::views::zip(hash_values, bin_indices))
        ibf.emplace(hash, seqan::hibf::bin_index{bin});

    auto agent = ibf.counting_agent();
    for (auto _ : state)
    {
        [[maybe_unused]] auto & res = agent.bulk_count(hash_values);
        benchmark::ClobberMemory();
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

BENCHMARK(emplace_benchmark)->Apply(arguments);
BENCHMARK(clear_benchmark)->Apply(arguments);
BENCHMARK(clear_range_benchmark)->Apply(arguments);
BENCHMARK(bulk_contains_benchmark)->Apply(arguments);
BENCHMARK(bulk_count_benchmark)->Apply(arguments);

BENCHMARK_MAIN();
