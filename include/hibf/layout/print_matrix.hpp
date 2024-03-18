// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <iostream>

#include <hibf/platform.hpp>

namespace seqan::hibf::layout
{

/*!\brief Helper function to print a matrix when debugging.
 * \ingroup hibf_layout
 */
template <typename matrix_type, typename matrix_value_type>
void print_matrix(matrix_type const & matrix,
                  size_t const row_bound,
                  size_t const column_bound,
                  matrix_value_type const inf)
{
    for (size_t i = 0; i < row_bound; ++i)
    {
        for (size_t j = 0; j < column_bound; ++j)
        {
            if (matrix[i][j] == inf)
                std::cerr << "inf\t";
            else
                std::cerr << matrix[i][j] << '\t';
        }
        std::cerr << '\n';
    }
    std::cerr << '\n';
}

} // namespace seqan::hibf::layout
