// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h> // for State, DoNotOptimize, Benchmark, BENCHMARK_DEFINE_F, BENCHMARK_...

#include <algorithm>  // for all_of, copy, __all_of_fn
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for identity
#include <memory>     // for allocator

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter

#include <sdsl/int_vector.hpp> // for operator-

using bitvector_t = seqan::hibf::interleaved_bloom_filter::binning_bitvector;

static void arguments(benchmark::internal::Benchmark * b)
{
    // 1,024 bits (1 KiB)
    b->Args({1024, 1024}); // No bit set
    b->Args({1024, 0});    // First bit set
#if 0
    b->Args({1024, 64});   // First bit in second word set
    b->Args({1024, 512});  // First bit set halfway

    // 8,192 bits (8 KiB)
    b->Args({8192, 8192}); // No bit set
    b->Args({8192, 0});    // First bit set
    b->Args({8192, 64});   // First bit in second word set
    b->Args({8192, 4096}); // First bit set halfway

    // 1,048,576 bits (1 MiB)
    b->Args({1LL << 20, 1LL << 20}); // No bit set
    b->Args({1LL << 20, 0});         // First bit set
    b->Args({1LL << 20, 64});        // First bit in second word set
    b->Args({1LL << 20, 1LL << 19}); // First bit set halfway
#endif
}

class all_zero : public benchmark::Fixture
{
public:
    void SetUp(benchmark::State const & state)
    {
        size_t const size_in_bits = state.range(0);
        size_t const first_set_bit = state.range(1);

        bitvector = bitvector_t(size_in_bits);
        if (first_set_bit < size_in_bits)
            bitvector[first_set_bit] = true;
    }

    void TearDown(benchmark::State const &)
    {
        bitvector.raw_data().clear();
    }

    bitvector_t const & get_bitvector() noexcept
    {
        return bitvector;
    }

private:
    bitvector_t bitvector{};
};

BENCHMARK_DEFINE_F(all_zero, std_all_of)(benchmark::State & state)
{
    bitvector_t const & bitvector = get_bitvector();

    for (auto _ : state)
    {
        bool result = std::all_of(bitvector.begin(),
                                  bitvector.end(),
                                  [](bool const value)
                                  {
                                      return !value;
                                  });
        benchmark::DoNotOptimize(result);
    }
}

BENCHMARK_DEFINE_F(all_zero, std_ranges_all_of)(benchmark::State & state)
{
    bitvector_t const & bitvector = get_bitvector();

    for (auto _ : state)
    {
        bool result = std::ranges::all_of(bitvector,
                                          [](bool const value)
                                          {
                                              return !value;
                                          });
        benchmark::DoNotOptimize(result);
    }
}

bool no_early_termination(bitvector_t const & bitvector) noexcept
{
    uint64_t const * const ptr = bitvector.raw_data().data();
    size_t const number_of_words{(bitvector.size() + 63u) >> 6};
    bool result{false};

    for (size_t i{}; i < number_of_words; ++i)
        result |= ptr[i];

    return !result;
}

BENCHMARK_DEFINE_F(all_zero, ptr_no_early_termination)(benchmark::State & state)
{
    bitvector_t const & bitvector = get_bitvector();

    for (auto _ : state)
    {
        bool result = no_early_termination(bitvector);
        benchmark::DoNotOptimize(result);
    }
}

bool with_early_termination(bitvector_t const & bitvector) noexcept
{
    uint64_t const * const ptr = bitvector.raw_data().data();
    size_t const number_of_words{(bitvector.size() + 63u) >> 6};
    bool result{false};

    for (size_t i{}; !result && i < number_of_words; ++i)
        result |= ptr[i];

    return !result;
}

BENCHMARK_DEFINE_F(all_zero, ptr_with_early_termination)(benchmark::State & state)
{
    bitvector_t const & bitvector = get_bitvector();

    for (auto _ : state)
    {
        bool result = with_early_termination(bitvector);
        benchmark::DoNotOptimize(result);
    }
}

BENCHMARK_REGISTER_F(all_zero, std_all_of)->Apply(arguments);
BENCHMARK_REGISTER_F(all_zero, std_ranges_all_of)->Apply(arguments);
BENCHMARK_REGISTER_F(all_zero, ptr_no_early_termination)->Apply(arguments);
BENCHMARK_REGISTER_F(all_zero, ptr_with_early_termination)->Apply(arguments);

BENCHMARK_MAIN();
