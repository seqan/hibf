// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-FileCopyrightText: 2013 Hideaki Ohno <hide.o.j55{at}gmail.com>
// SPDX-License-Identifier: BSD-3-Clause AND MIT

#include <algorithm>  // for max, __count_fn, __fill_fn, count, fill
#include <array>      // for array
#include <bit>        // for countl_zero
#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t, uint8_t, uint32_t
#include <cmath>      // for log
#include <cstddef>    // for size_t
#include <functional> // for identity
#include <iostream>   // for basic_ostream, basic_istream, basic_istream::read, basic_ostre...
#include <stdexcept>  // for runtime_error, invalid_argument
#include <utility>    // for addressof, swap
#include <vector>     // for vector

#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/sketch/hyperloglog.hpp>        // for hyperloglog

#include <simde/x86/avx.h>  // for simde__m256i
#include <simde/x86/avx2.h> // for simde_mm256_max_epu8

namespace seqan::hibf::sketch
{

hyperloglog::hyperloglog(uint8_t const bits) : bits{bits}, size{1ULL << bits}, data(size, 0u)
{
    if (bits < 5u || bits > 32u)
        throw std::invalid_argument("[HyperLogLog] bit width must be in the range [5,32].");

    double correction_factor{}; // bias-correction, "alpha"

    switch (size)
    {
    case 32:
        correction_factor = 0.697;
        break;
    case 64:
        correction_factor = 0.709;
        break;
    default:
        correction_factor = 0.7213 / (1.0 + 1.079 / size);
        break;
    }

    normalization_factor = correction_factor * size * size;
    // Last `bits` bits are 1.
    rank_mask = (1ULL << bits) - 1u;
}

// See https://github.com/wangyi-fudan/wyhash
// Simpler than murmur3: https://godbolt.org/z/bz7fd4aYz
[[nodiscard]] inline uint64_t wyhash(uint64_t const value) noexcept
{
    __uint128_t result = value;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}

void hyperloglog::add(uint64_t const value)
{
    uint64_t const hash = wyhash(value);
    // The first bits bits are used to distribute the leading zero counts over data.
    uint64_t const index = hash >> (64 - bits);
    // rank_mask ensures that the lzcount is at most 64 - bits.
    uint8_t const rank = std::countl_zero((hash << bits) | rank_mask) + 1;
    data[index] = std::max(rank, data[index]);
}

double hyperloglog::estimate() const
{
    float sum = 0.0;

#pragma omp simd
    for (size_t i = 0; i < size; ++i)
        sum += expectation_values[data[i]];

    double estimate = normalization_factor / sum;

    // Small value correction: linear counting of zeros
    if (estimate <= 2.5 * size)
    {
        uint32_t const zeros = std::ranges::count(data, uint8_t{});
        if (zeros != 0u)
            estimate = size * std::log(static_cast<double>(size) / zeros);
    }

    return estimate;
}

void hyperloglog::merge(hyperloglog const & other)
{
    assert(size == other.size);

    simde__m256i * const it = reinterpret_cast<simde__m256i *>(data.data());
    simde__m256i const * const other_it = reinterpret_cast<simde__m256i const *>(other.data.data());

    // We can do 256 bits = 32 bytes at once.
    // We store `uint8_t`, so `size` is the size in bytes.
    // Hence, we need to do `size / 32` iterations.
    for (size_t i = 0; i < size / 32; ++i)
    {
        it[i] = simde_mm256_max_epu8(it[i], other_it[i]);
    }
}

double hyperloglog::merge_and_estimate(hyperloglog const & other)
{
    merge(other);
    return estimate();
}

void hyperloglog::reset()
{
    std::ranges::fill(data, 0u);
}

void hyperloglog::store(std::ostream & os) const
{
    assert(data.size() == size);

    char const * const bits_ptr = reinterpret_cast<char const *>(std::addressof(bits));
    os.write(bits_ptr, sizeof(bits));

    char const * const data_ptr = reinterpret_cast<char const *>(data.data());
    os.write(data_ptr, sizeof(data[0]) * size);

    os.flush();
    if (os.fail())
    {
        throw std::runtime_error("[HyperLogLog] Failed to store a HyperLogLog sketch to a file.");
    }
}

void hyperloglog::load(std::istream & is)
{
    try
    {
        uint8_t restore_bits{};
        char * const bits_ptr = reinterpret_cast<char *>(std::addressof(restore_bits));
        is.read(bits_ptr, sizeof(restore_bits));

        hyperloglog restore_hll{restore_bits}; // Constructor might throw std::invalid_argument

        char * const data_ptr = reinterpret_cast<char *>(restore_hll.data.data());
        is.read(data_ptr, sizeof(data[0]) * restore_hll.size);

        if (is.fail())
        {
            throw std::runtime_error("[HyperLogLog] Failed to load a HyperLogLog sketch from a file: I/O error.");
        }
        std::swap(*this, restore_hll);
    }
    catch (std::invalid_argument const & err)
    {
        throw std::runtime_error("[HyperLogLog] Failed to load a HyperLogLog sketch from a file: Invalid bit_width.");
    }
}

} // namespace seqan::hibf::sketch
