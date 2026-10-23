// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-FileCopyrightText: 2013 Hideaki Ohno <hide.o.j55{at}gmail.com>
// SPDX-License-Identifier: BSD-3-Clause AND MIT

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

[[nodiscard]] uint64_t wyhash(uint64_t const value) noexcept;

} // namespace seqan::hibf::sketch
