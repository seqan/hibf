// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for AddCustomContext, State, Benchmark, Counter, DoNotOptimize, BEN...

#include <cmath>      // for log, ceil, exp
#include <cstddef>    // for size_t
#include <functional> // for function
#include <ranges>     // for iota_view, all_t, view_interface, __fn, iota, views
#include <string>     // for to_string, basic_string

#include <hibf/config.hpp>                   // for insert_iterator, config
#include <hibf/contrib/std/chunk_view.hpp>   // for chunk_view, chunk, chunk_fn
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_count, bin_index, bin_size, hash_...
#include <hibf/misc/divide_and_ceil.hpp>     // for divide_and_ceil
#include <hibf/platform.hpp>                 // for HIBF_HAS_AVX512
#include <hibf/test/bytes.hpp>               // for operator""_MiB

using namespace seqan::hibf::test::literals;
static constexpr size_t total_ibf_size_in_bytes{1_MiB};
static constexpr size_t number_of_hash_functions{2u};
static constexpr double false_positive_rate{0.05};

inline benchmark::Counter ibf_size(size_t const bit_size)
{
    return benchmark::Counter(bit_size / 8, benchmark::Counter::kDefaults, benchmark::Counter::OneK::kIs1024);
}

// This computes how many elements need to be inserted into the IBF to achieve the desired false positive rate for the
// given size.
// The `number_of_elements` many generated values are used for both constructing and querying the IBF.
static /* cmath not constexpr in libc++ */ size_t number_of_elements = []()
{
    size_t const bits = 8u * total_ibf_size_in_bytes;
    double const numerator = -std::log(1 - std::exp(std::log(false_positive_rate) / number_of_hash_functions)) * bits;
    return std::ceil(numerator / number_of_hash_functions);
}();

static auto get_value(size_t const bins)
{
    size_t const chunk_size = seqan::hibf::divide_and_ceil(number_of_elements, bins);
    return seqan::stl::views::chunk(std::views::iota(size_t{}, number_of_elements), chunk_size);
}

void manual_construct(::benchmark::State & state)
{
    size_t const bins = state.range(0);
    size_t const bits = 8u * total_ibf_size_in_bytes / bins;

    auto values = get_value(bins);

    for (auto _ : state)
    {
        seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{bins},
                                                  seqan::hibf::bin_size{bits},
                                                  seqan::hibf::hash_function_count{number_of_hash_functions}};

        for (size_t bin_index = 0u; bin_index < bins; ++bin_index)
        {
            for (auto value : values[bin_index])
                ibf.emplace(value, seqan::hibf::bin_index{bin_index});
        }

        state.counters["IBF_size"] = ibf_size(ibf.bit_size());

        benchmark::DoNotOptimize(ibf);
    }
}

void config_construct(::benchmark::State & state)
{
    size_t const bins = state.range(0);

    auto values = get_value(bins);

    seqan::hibf::config config{.input_fn =
                                   [&values](size_t const user_bin_id, seqan::hibf::insert_iterator && it)
                               {
                                   for (auto const value : values[user_bin_id])
                                       it = value;
                               },
                               .number_of_user_bins = bins,
                               .number_of_hash_functions = number_of_hash_functions,
                               .maximum_fpr = false_positive_rate};

    for (auto _ : state)
    {
        seqan::hibf::interleaved_bloom_filter ibf{config};

        state.counters["IBF_size"] = ibf_size(ibf.bit_size());

        benchmark::DoNotOptimize(ibf);
    }
}

void config_and_max_construct(::benchmark::State & state)
{
    size_t const bins = state.range(0);

    auto values = get_value(bins);
    size_t const max_bin_size = values[0].size();

    seqan::hibf::config config{.input_fn =
                                   [&values](size_t const user_bin_id, seqan::hibf::insert_iterator && it)
                               {
                                   for (auto const value : values[user_bin_id])
                                       it = value;
                               },
                               .number_of_user_bins = bins,
                               .number_of_hash_functions = number_of_hash_functions,
                               .maximum_fpr = false_positive_rate};

    for (auto _ : state)
    {
        seqan::hibf::interleaved_bloom_filter ibf{config, max_bin_size};

        state.counters["IBF_size"] = ibf_size(ibf.bit_size());

        benchmark::DoNotOptimize(ibf);
    }
}

BENCHMARK(manual_construct)->RangeMultiplier(2)->Range(64, 1024);
BENCHMARK(config_construct)->RangeMultiplier(2)->Range(64, 1024);
BENCHMARK(config_and_max_construct)->RangeMultiplier(2)->Range(64, 1024);

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
