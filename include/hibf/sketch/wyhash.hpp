// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::sketch::wyhash.
 * See https://github.com/wangyi-fudan/wyhash
 * Simpler than murmur3: https://godbolt.org/z/bz7fd4aYz
 */

#pragma once

#include <cstdint> // for uint64_t

#include <hibf/platform.hpp>

namespace seqan::hibf::sketch
{

/*!\brief Mixes the bits of `value`: the 128-bit product with a 64-bit constant, folded to 64 bits via XOR.
 * \ingroup hibf_sketch
 */
[[nodiscard]] constexpr uint64_t wyhash(uint64_t const value) noexcept
{
    __uint128_t result = value;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}

} // namespace seqan::hibf::sketch
