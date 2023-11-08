// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h> // for Benchmark, State, BENCHMARK_CAPTURE, DoNotOptimize, BENCHMARK_MAIN

#include <algorithm> // for __generate_fn, generate
#include <cinttypes> // for int32_t, uint8_t
#include <cstddef>   // for size_t
#include <filesystem>
#include <fstream>
#include <memory>      // for allocator
#include <random>      // for uniform_int_distribution, mt19937_64
#include <type_traits> // for invoke_result_t
#include <utility>     // for move, pair

#include <hibf/misc/bit_vector.hpp>     // for bit_vector
#include <hibf/test/sandboxed_path.hpp> // for sandboxed_path, operator/
#include <hibf/test/tmp_directory.hpp>  // for tmp_directory

#include <cereal/archives/binary.hpp> // for BinaryInputArchive, BinaryOutputArchive

seqan::hibf::bit_vector generate_bit_vector(size_t const size_in_bits, size_t const seed = 0u)
{
    std::mt19937_64 engine{seed};
    std::uniform_int_distribution<uint8_t> dist{0u, 1u};

    auto gen = [&dist, &engine]()
    {
        return dist(engine);
    };
    seqan::hibf::bit_vector vec(size_in_bits);
    std::ranges::generate(vec, gen);

    return vec;
}

static seqan::hibf::test::tmp_directory tmp{};

void load_bit_vector(benchmark::State & state)
{
    size_t const size_in_bits = 1ULL << state.range(0);
    auto const filename = tmp.path() / std::to_string(state.range(0));

    if (!std::filesystem::exists(filename))
    {
        seqan::hibf::bit_vector const vector = generate_bit_vector(size_in_bits);

        std::ofstream output_stream{filename, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{output_stream};
        oarchive(vector);
    }

    // Substract 8 bytes for serialised size_ member.
    size_t const filesize_in_bytes = std::filesystem::file_size(filename) - 8u;

    if (size_in_bits / 8u != filesize_in_bytes)
        throw std::logic_error{"Actual and expected file size differ."};

    for (auto _ : state)
    {
        seqan::hibf::bit_vector vector{};
        {
            std::ifstream input_stream{filename, std::ios::binary};
            cereal::BinaryInputArchive iarchive{input_stream};
            iarchive(vector);
        }
        benchmark::DoNotOptimize(vector);
    }

    state.SetBytesProcessed(static_cast<int64_t>(state.iterations()) * filesize_in_bytes);
    // Use `filesize_in_bytes` to make benchmark show `1Mi` instead of `1.00001Mi`.
    state.counters["filesize"] =
        benchmark::Counter(filesize_in_bytes, benchmark::Counter::kDefaults, benchmark::Counter::OneK::kIs1024);
}

// Benchmark ranges are `int32_t`, i.e to small to represent 4 GiB or more.
// Hence we use `2^x`.
// 13 -> 1KiB
// 23 -> 1MiB
// 33 -> 1GiB
// We use a small number here for the unit tests.
// Sizes of 10s of GiB would be more interesting for actual benchmarking.
static constexpr int32_t min_range = 13;
static constexpr int32_t max_range = 13;

BENCHMARK(load_bit_vector)->Range(min_range, max_range);

BENCHMARK_MAIN();
