// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cassert>    // for assert
#include <filesystem> // for remove
#include <fstream>    // for char_traits, basic_ofstream, basic_ostream, operator<<, ofstream

#include <hibf/test/sandboxed_path.hpp> // for operator/, sandboxed_path
#include <hibf/test/tmp_directory.hpp>  // for tmp_directory

int main()
{
    // create a directory folder
    seqan::hibf::test::tmp_directory tmp{};

    // Some function that should creates temporary files and removes them again
    {
        std::ofstream ofs{tmp.path() / "somefile.txt"};
        ofs << "Hello World!";
        ofs.close();

        std::filesystem::remove(tmp.path() / "somefile.txt");
    }

    // check that everything was cleaned up properly
    assert(tmp.empty());
}
