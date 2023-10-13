// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cinttypes> // for int64_t, int16_t, int32_t, int8_t, uint16_t, uint32_t, uint64_t
#include <iostream>  // for cout, ostream
#include <vector>    // for vector

#include <hibf/interleaved_bloom_filter.hpp> // for counting_vector, binning_bitvector

namespace seqan::hibf
{

// Why a function object and not free `functions void print(...)` ?
// Disables ADL and requires fully qualified name outside the seqan::hibf namespace.
// With a free function, both seqan::hibf::print and print (found via ADL) would work.
// Print is a common name and a free function might cause clashes.
// A free function for `std::vector<int64_t>` might also cause problems.
struct print_t
{
    void operator()(seqan::hibf::binning_bitvector const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<uint8_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<uint16_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<uint32_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<uint64_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<int8_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<int16_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<int32_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(seqan::hibf::counting_vector<int64_t> const & vector, std::ostream & stream = std::cout) const;
    void operator()(std::vector<int64_t> const & vector, std::ostream & stream = std::cout) const;
};

static inline constexpr auto print = print_t{};

} // namespace seqan::hibf
