// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan::hibf::md_vector.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <vector> // for vector

#include <hibf/platform.hpp>

namespace seqan::hibf
{

template <typename value_t>
struct md_vector : public std::vector<std::vector<value_t>>
{
    using base_t = std::vector<std::vector<value_t>>;
    using base_t::base_t;
    using base_t::operator[];
#if defined(__cpp_explicit_this_parameter) && __cpp_explicit_this_parameter >= 202110L
    decltype(auto) operator[](this auto & self, size_t const x, size_t const y)
    {
        return self[x][y];
    }
#else
    value_t & operator[](size_t const x, size_t const y)
    {
        return (*this)[x][y];
    }
    value_t const & operator[](size_t const x, size_t const y) const
    {
        return (*this)[x][y];
    }
#endif
};

} // namespace seqan::hibf
