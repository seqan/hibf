// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <hibf/misc/unreachable.hpp> // for unreachable

int foo(int const i)
{
    // The compiler will not generate the default branch.
    // Note that an input of any `i` other than `0` and `1` is undefined behavior!
    switch (i)
    {
    case 0:
        return -5;
    case 1:
        return 3;
    default:
        seqan::hibf::unreachable();
    }
}
