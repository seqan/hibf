// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <cassert>     // for assert
#include <cmath>       // for ceil, sqrt
#include <functional>  // for function
#include <iostream>    // for operator<<, basic_ostream, basic_istream, getline, stringstream
#include <sstream>     // for basic_stringstream
#include <stdexcept>   // for invalid_argument
#include <string>      // for char_traits, string
#include <string_view> // for operator==, basic_string_view, string_view

#include <cereal/archives/json.hpp> // for JSONInputArchive, JSONOutputArchive
#include <cereal/cereal.hpp>        // for make_nvp, InputArchive, OutputArchive

#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/misc/insert_iterator.hpp> // for insert_iterator
#include <hibf/sketch/hyperloglog.hpp>

namespace seqan::hibf
{

insert_iterator & insert_iterator::operator=(uint64_t const value) noexcept
{
    assert(ptr != nullptr);

    switch (type)
    {
    case data_type::unordered_set:
        static_cast<set_t *>(ptr)->emplace(value);
        break;
    case data_type::sketch:
        static_cast<sketch_t *>(ptr)->add(value);
        break;
    case data_type::ibf:
        static_cast<ibf_t *>(ptr)->emplace(value, static_cast<bin_index>(ibf_bin_index));
        break;
    default:
#ifndef NDEBUG
        assert(false);
#else
        __builtin_unreachable();
#endif
        // assert(type == data_type::function);
        // static_cast<function_t *>(ptr)->operator()(value);
    }
    return *this;
}

} // namespace seqan::hibf
