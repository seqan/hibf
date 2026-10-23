// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstdint> // for uint64_t

#include <hibf/sketch/wyhash.hpp> // for wyhash

namespace seqan::hibf::sketch
{

// See https://github.com/wangyi-fudan/wyhash
// Simpler than murmur3: https://godbolt.org/z/bz7fd4aYz
[[nodiscard]] uint64_t wyhash(uint64_t const value) noexcept
{
    __uint128_t result = value;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}

} // namespace seqan::hibf::sketch
