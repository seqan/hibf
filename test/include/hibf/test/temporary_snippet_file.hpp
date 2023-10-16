// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
    temporary_snippet_file(std::string_view const & file_name, content_t &&... content) :
        path_{tmp_folder.path() / file_name}
    {
        std::filesystem::current_path(tmp_folder.path());

        if constexpr (sizeof...(content_t) > 0)
        {
            std::ofstream file{path_};
            (file << ... << std::forward<content_t>(content));
        }
    }

    std::filesystem::path path() const
    {
        return path_;
    }

private:
    static inline seqan::hibf::test::tmp_directory const tmp_folder{};
    std::filesystem::path path_{};
};

} // namespace seqan::hibf::test
