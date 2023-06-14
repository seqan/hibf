// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for AssertionResult, TestInfo, Message, TEST_F, TestPartResult, EXPECT_EQ

#include <filesystem> // for remove, path, temp_directory_path, permissions, perm_options, perms
#include <iosfwd>     // for fstream, ios, ofstream
#include <memory>     // for allocator

#include <hibf/test/file_access.hpp> // for write_access, read_access

struct file_access_test : public ::testing::Test
{
    static std::filesystem::path create_file(char const * const file_name)
    {
        std::filesystem::path file = std::filesystem::temp_directory_path();
        file /= file_name;

        std::ofstream str{file};

        return file;
    }

    static std::filesystem::path create_directory(char const * const directory_name)
    {
        std::filesystem::path directory = std::filesystem::temp_directory_path();
        directory /= directory_name;

        std::filesystem::create_directory(directory);

        return directory;
    }

    static void remove_read_permission(std::filesystem::path const & path)
    {
        std::filesystem::permissions(path, std::filesystem::perms::owner_read, std::filesystem::perm_options::remove);
    }

    static void remove_write_permission(std::filesystem::path const & path)
    {
        std::filesystem::permissions(path, std::filesystem::perms::owner_write, std::filesystem::perm_options::remove);
    }

    // If we have root permissions, we can write to the file even when we do not have write permissions.
    static bool is_root()
    {
        std::filesystem::path tmp_file = std::filesystem::temp_directory_path();
        tmp_file /= "hibf_test_permissions_is_root.txt";

        {
            std::fstream stream{tmp_file, std::ios::out};
        }

        remove_write_permission(tmp_file);

        bool result{};

        {
            std::fstream stream;
            stream.open(tmp_file, std::ios::out);
            result = !stream.fail();
        }

        std::filesystem::remove(tmp_file);

        return result;
    }
};

TEST_F(file_access_test, file_read_access_granted)
{
    auto path = create_file("hibf_test_permissions_file_read_access_granted");
    EXPECT_TRUE(hibf::test::read_access(path));
    std::filesystem::remove(path);
}

TEST_F(file_access_test, file_read_access_revoked)
{
    auto path = create_file("hibf_test_permissions_file_read_access_revoked");
    remove_read_permission(path);
    EXPECT_EQ(hibf::test::read_access(path), is_root());
    std::filesystem::remove(path);
}

TEST_F(file_access_test, file_write_access_granted)
{
    auto path = create_file("hibf_test_permissions_file_write_access_granted");
    EXPECT_TRUE(hibf::test::write_access(path));
    std::filesystem::remove(path);
}

TEST_F(file_access_test, file_write_access_revoked)
{
    auto path = create_file("hibf_test_permissions_file_write_access_revoked");
    remove_write_permission(path);
    EXPECT_EQ(hibf::test::write_access(path), is_root());
    std::filesystem::remove(path);
}

TEST_F(file_access_test, directory_write_access_granted)
{
    auto path = create_directory("hibf_test_permissions_directory_write_access_granted");
    EXPECT_TRUE(hibf::test::write_access(path));
    std::filesystem::remove(path);
}

TEST_F(file_access_test, directory_write_access_revoked)
{
    auto path = create_directory("hibf_test_permissions_directory_write_access_revoked");
    remove_write_permission(path);
    EXPECT_EQ(hibf::test::write_access(path), is_root());
    std::filesystem::remove(path);
}
