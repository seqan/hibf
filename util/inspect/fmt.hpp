// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <fmt/color.h>

#include <sharg/detail/terminal.hpp>

namespace seqan::hibf::util
{

struct styles
{
    static inline bool const cout_is_terminal = sharg::detail::stdout_is_terminal();
    static inline bool const cerr_is_terminal = sharg::detail::stderr_is_terminal();

    static bool use_style(FILE * stream)
    {
        return (stream == stdout && cout_is_terminal) || (stream == stderr && cerr_is_terminal);
    }

    static fmt::text_style color(fmt::color const color, FILE * stream = stdout)
    {
        if (use_style(stream))
            return fmt::fg(color);
        else
            return fmt::text_style{};
    }

    static fmt::text_style emphasis(fmt::emphasis const emphasis, FILE * stream = stdout)
    {
        if (use_style(stream))
            return emphasis;
        else
            return fmt::text_style{};
    }
};

} // namespace seqan::hibf::util
