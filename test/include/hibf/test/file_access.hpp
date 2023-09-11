// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

/*!\file
 * \brief Detects read and write access for a path.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem> // for path, is_directory, remove
#include <fstream>    // for fstream, ios

#include <hibf/platform.hpp>

namespace seqan::hibf::test
{

//!\brief Checks wheter there is read access to a path.
inline bool read_access(std::filesystem::path const & file)
{
    std::fstream stream;
    stream.open(file, std::ios::in);
    return !stream.fail();
}

//!\brief Checks wheter there is write access to a path.
inline bool write_access(std::filesystem::path const & file)
{
    if (std::filesystem::is_directory(file))
    {
        std::filesystem::path test_file{file};
        test_file /= "hibf_test_write_access";
        std::fstream stream;
        stream.open(test_file, std::ios::out);
        bool result = !stream.fail();
        if (result)
        {
            stream.close();
            std::filesystem::remove(test_file);
        }
        return result;
    }
    else
    {
        std::fstream stream;
        stream.open(file, std::ios::out);
        return !stream.fail();
    }
}

} // namespace seqan::hibf::test
