// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <iostream>

#include <hibf/misc/print.hpp>

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

void print_t::operator()(seqan::hibf::binning_bitvector const & vector, std::ostream & stream) const
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

void print_t::operator()(std::vector<int64_t> const & vector, std::ostream & stream) const
{
    print_impl(vector, stream);
}

} // namespace seqan::hibf
