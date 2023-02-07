// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <library_template/test/temporary_snippet_file.hpp>

TEST(temporary_snippet_file, no_content)
{
    std::filesystem::path const original_path = std::filesystem::current_path();
    library_template::test::temporary_snippet_file tmp_file{"test.txt"};
    std::filesystem::path const new_path = std::filesystem::current_path();

    EXPECT_FALSE(std::filesystem::exists(new_path / "test.txt"));
    EXPECT_NE(original_path, new_path); // Only check this once, because temporary_snippet_file sets current_path.
}

TEST(temporary_snippet_file, with_content)
{
    std::filesystem::path const path = std::filesystem::current_path();
    library_template::test::temporary_snippet_file tmp_file{"test.txt", "some content\n", "more"};

    EXPECT_TRUE(std::filesystem::exists(path / "test.txt"));

    std::ifstream file{path / "test.txt"};
    std::stringstream file_buffer;
    file_buffer << file.rdbuf();
    EXPECT_EQ(file_buffer.str(), "some content\nmore");
}
