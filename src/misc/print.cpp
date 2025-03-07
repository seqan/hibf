// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cstddef>     // for size_t
#include <cstdint>     // for uint16_t, int16_t, int8_t, uint64_t, uint8_t, int32_t, int64_t, uin...
#include <iostream>    // for ostream, basic_ostream, operator<<
#include <ranges>      // for range_value_t, empty
#include <type_traits> // for conditional
#include <vector>      // for vector

#include <hibf/misc/bit_vector.hpp>      // for bit_vector
#include <hibf/misc/counting_vector.hpp> // for counting_vector
#include <hibf/misc/print.hpp>           // for print_t

namespace seqan::hibf
{

template <typename value_t>
struct printable_value_t
{
    using value = value_t;
};

template <>
struct printable_value_t<uint8_t>
{
    using value = uint16_t;
};

template <>
struct printable_value_t<int8_t>
{
    using value = int16_t;
};

template <>
struct printable_value_t<bool>
{
    using value = uint16_t;
};

template <typename vector_t>
void print_impl(vector_t const & vector, std::ostream & stream)
{
    using value_t = printable_value_t<std::ranges::range_value_t<vector_t>>::value;

    stream << '[';

    if (!std::ranges::empty(vector))
    {
        for (size_t i = 0u; i < vector.size() - 1u; ++i)
            stream << static_cast<value_t>(vector[i]) << ',';
        stream << static_cast<value_t>(vector[vector.size() - 1u]);
    }

    stream << "]\n";
}

void print_t::operator()(seqan::hibf::bit_vector const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<uint8_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<uint16_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<uint32_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<uint64_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<int8_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<int16_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<int32_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(seqan::hibf::counting_vector<int64_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

void print_t::operator()(std::vector<uint64_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

} // namespace seqan::hibf
