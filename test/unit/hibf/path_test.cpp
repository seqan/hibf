// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h> // for Test, TestInfo, TEST

#include <filesystem> // for path
#include <utility>    // for move

#include <cereal/specialize.hpp> // for specialization

#include <hibf/cereal/path.hpp> // IWYU pragma: keep
#include <hibf/test/cereal.hpp> // for test_serialisation

TEST(path_test, serialisation)
{
    std::filesystem::path path{"/some/random/path.txt"};
    seqan::hibf::test::test_serialisation(std::move(path));
}
