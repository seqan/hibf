// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem> // for path
#include <string>     // for string

#include <cereal/macros.hpp>       // for CEREAL_LOAD_FUNCTION_NAME, CEREAL_SAVE_FUNCTION_NAME
#include <cereal/types/string.hpp> // IWYU pragma: keep

#include <hibf/platform.hpp>

namespace cereal
{

template <typename archive_t>
void CEREAL_SAVE_FUNCTION_NAME(archive_t & archive, std::filesystem::path const & path)
{
    std::string const str{path.string()};
    archive(str);
}

template <typename archive_t>
void CEREAL_LOAD_FUNCTION_NAME(archive_t & archive, std::filesystem::path & path)
{
    std::string str;
    archive(str);
    path.assign(str);
}

} // namespace cereal
