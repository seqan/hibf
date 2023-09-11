// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

/*!\file
 * \brief Utility for creating snippet files.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>  // for current_path
#include <fstream>     // for ofstream
#include <string_view> // for string_view
#include <utility>     // for forward

#include <hibf/test/sandboxed_path.hpp> // for sandboxed_path, operator/
#include <hibf/test/tmp_directory.hpp>  // for tmp_directory

namespace seqan::hibf::test
{

class temporary_snippet_file
{
public:
    temporary_snippet_file() = default;
    temporary_snippet_file(temporary_snippet_file const &) = default;
    temporary_snippet_file(temporary_snippet_file &&) = default;
    temporary_snippet_file & operator=(temporary_snippet_file const &) = default;
    temporary_snippet_file & operator=(temporary_snippet_file &&) = default;
    ~temporary_snippet_file() = default;

    template <typename... content_t>
    temporary_snippet_file(std::string_view const & file_name, content_t &&... content)
    {
        std::filesystem::current_path(tmp_folder.path());

        if constexpr (sizeof...(content_t) > 0)
        {
            std::ofstream file{tmp_folder.path() / file_name};
            (file << ... << std::forward<content_t>(content));
        }
    }

private:
    static inline seqan::hibf::test::tmp_directory const tmp_folder{};
};

} // namespace seqan::hibf::test
