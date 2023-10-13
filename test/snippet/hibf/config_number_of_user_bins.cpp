// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cstddef>    // for size_t
#include <functional> // for function

#include <hibf/config.hpp> // for config, insert_iterator

int main()
{
    auto my_input = [&](size_t const /* user_bin_id */, seqan::hibf::insert_iterator it) // fixed parameters!
    {
        it = 42; // assign something that is convertible to uint64_t
    };

    seqan::hibf::config config{.input_fn = my_input, .number_of_user_bins = 12};
}
