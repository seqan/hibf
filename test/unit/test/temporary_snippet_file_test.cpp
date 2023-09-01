// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, AssertionResult, Message, TestInfo, TestPartResult

#include <filesystem> // for current_path, operator/, path, exists
#include <fstream>    // for char_traits, basic_ostream::operator<<, basic_ifstream, basi...
#include <memory>     // for allocator
#include <sstream>    // for basic_stringstream

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
    std::filesystem::path const path = std::filesystem::current_path();
    seqan::hibf::test::temporary_snippet_file tmp_file{"test.txt", "some content\n", "more"};
    std::filesystem::path const new_path = std::filesystem::current_path();

    EXPECT_TRUE(std::filesystem::exists(new_path / "test.txt"));

    std::ifstream file{new_path / "test.txt"};
    std::stringstream file_buffer;
    file_buffer << file.rdbuf();
    EXPECT_EQ(file_buffer.str(), "some content\nmore");
}
