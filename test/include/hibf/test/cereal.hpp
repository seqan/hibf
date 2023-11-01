// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides cereal functionality for tests.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h> // for AssertionResult, EXPECT_TRUE, Message, TestPartResult

#include <fstream>     // for basic_ios, ios, basic_ifstream, basic_ofstream, ifstream, ofs...
#include <type_traits> // for remove_cvref_t

#include <hibf/cereal/concepts.hpp>     // for cereal_input_archive, cereal_output_archive
#include <hibf/test/sandboxed_path.hpp> // for sandboxed_path, operator/
#include <hibf/test/tmp_directory.hpp>  // for tmp_directory

#include <cereal/archives/binary.hpp>          // for BinaryInputArchive, BinaryOutputArchive
#include <cereal/archives/json.hpp>            // for JSONInputArchive, JSONOutputArchive
#include <cereal/archives/portable_binary.hpp> // for PortableBinaryInputArchive, PortableBinaryOutputArchive
#include <cereal/archives/xml.hpp>             // for XMLInputArchive, XMLOutputArchive

namespace seqan::hibf::test
{

/*!\brief Tests if an object is serialiseable for a specific cereal archive type.
 * \tparam in_archive_t  Type of the cereal input archive. Must model seqan::hibf::cereal_input_archive.
 * \tparam out_archive_t Type of the cereal output archive. Must model seqan::hibf::cereal_output_archive.
 * \tparam value_t       The type to cerealise.
 * \param value The object to cerealise.
 */
template <cereal_input_archive in_archive_t, cereal_output_archive out_archive_t, typename value_t>
void test_serialisation(value_t && value)
{
    tmp_directory tmp{};
    auto const filename = tmp.path() / "cereal_test";

    {
        std::ofstream output_stream{filename, std::ios::binary};
        out_archive_t oarchive{output_stream};
        oarchive(value);
    }

    {
        std::remove_cvref_t<value_t> value_from_archive{};
        std::ifstream input_stream{filename, std::ios::binary};
        in_archive_t iarchive{input_stream};
        iarchive(value_from_archive);
        EXPECT_TRUE(value == value_from_archive);
    }
}

/*!\brief Tests if an object is serialiseable for all cereal archive types.
 * \tparam value_t The type to serialise.
 * \param value The object to serialise.
 */
template <typename value_t>
void test_serialisation(value_t && value)
{
    test_serialisation<cereal::BinaryInputArchive, cereal::BinaryOutputArchive>(value);
    test_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(value);
    test_serialisation<cereal::JSONInputArchive, cereal::JSONOutputArchive>(value);
    test_serialisation<cereal::XMLInputArchive, cereal::XMLOutputArchive>(value);
}

} // namespace seqan::hibf::test
