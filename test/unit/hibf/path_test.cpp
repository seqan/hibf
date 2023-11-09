// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, TestInfo, TEST

#include <filesystem> // for operator==, path
#include <memory>     // for allocator
#include <utility>    // for move

#include <hibf/test/cereal.hpp> // for test_serialisation

#include <cereal/specialize.hpp> // for specialization

TEST(path_test, serialisation)
{
    std::filesystem::path path{"/some/random/path.txt"};
    seqan::hibf::test::test_serialisation(std::move(path));
}
