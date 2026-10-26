// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::sketch::wyhash.
 */

#pragma once

#include <cstdint> // for uint64_t

#include <hibf/platform.hpp>

namespace seqan::hibf::sketch
{

/*!\brief Mixes the bits of `value`: the 128-bit product with a 64-bit constant, folded to 64 bits via XOR.
 * \ingroup hibf_sketch
 * \details
 * `wyhash(0) == 0`: the value 0 (e.g. the k-mer `AA...A` in a 2-bit encoding) always becomes the smallest possible
 * hash. It is therefore the smallest value of the first MinHash sketch of every user bin that contains it, and it sets
 * register 0 of a HyperLogLog sketch to the maximal rank. To avoid this, XOR `value` with a constant before the
 * multiplication; the input that is mapped to 0 is then this constant. This changes every hash, so sketches computed
 * before and after the change are not comparable.
 * \sa https://github.com/wangyi-fudan/wyhash
 * \sa Simpler than murmur3: https://godbolt.org/z/bz7fd4aYz
 */
[[nodiscard]] constexpr uint64_t wyhash(uint64_t const value) noexcept
{
    __uint128_t result = value;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}

} // namespace seqan::hibf::sketch
