// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, AssertionResult, Message, TestInfo, TestPartResult

#include <filesystem>  // for path, current_path, operator/, exists
#include <fstream>     // for char_traits, basic_ifstream, basic_filebuf, basic_ostream
#include <sstream>     // for basic_stringstream
#include <string>      // for allocator, basic_string
#include <string_view> // for basic_string_view

#include <hibf/test/temporary_snippet_file.hpp> // for temporary_snippet_file

TEST(temporary_snippet_file, no_content)
{
    std::filesystem::path const original_path = std::filesystem::current_path();
    seqan::hibf::test::temporary_snippet_file tmp_file{"test.txt"};
    std::filesystem::path const new_path = std::filesystem::current_path();

    EXPECT_FALSE(std::filesystem::exists(new_path / "test.txt"));
    EXPECT_NE(original_path, new_path); // Only check this once, because temporary_snippet_file sets current_path.
}

TEST(temporary_snippet_file, with_content)
{
    seqan::hibf::test::temporary_snippet_file tmp_file{"test.txt", "some content\n", "more"};
    std::filesystem::path const new_path = std::filesystem::current_path();

    EXPECT_TRUE(std::filesystem::exists(new_path / "test.txt"));

    std::ifstream file{new_path / "test.txt"};
    std::stringstream file_buffer;
    file_buffer << file.rdbuf();
    EXPECT_EQ(file_buffer.str(), "some content\nmore");
}
