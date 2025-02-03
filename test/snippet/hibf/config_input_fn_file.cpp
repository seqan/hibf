// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t
#include <filesystem>  // for path
#include <fstream>     // for basic_ifstream, getline, ifstream
#include <functional>  // for function
#include <string>      // for basic_string, string
#include <string_view> // for basic_string_view
#include <vector>      // for vector

#include <hibf/test/temporary_snippet_file.hpp> // for temporary_snippet_file

seqan::hibf::test::temporary_snippet_file file1{"file1.fa", "ACGT"};
seqan::hibf::test::temporary_snippet_file file2{"file2.fa", "ACGT"};

//![main]
#include <hibf/config.hpp> // for insert_iterator, config

int main()
{
    std::vector<std::filesystem::path> filenames{file1.path(), file2.path()};

    auto my_input = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        std::ifstream infile{filenames[user_bin_id]};
        std::string line;
        std::getline(infile, line); // assume there is a sequence in the first line
        // Look at https://docs.seqan.de/seqan3/3-master-user/group__io__sequence__file.html for e.g. FASTA File I/O

        for (size_t i = 0; i < line.size() - 1; ++i)
        {
            // compute 2-mer hashes based on the character value
            uint64_t hash = 4 * line[i] + line[i + 1];
            // You can also look at the seqan3::kmer_hash view for hashing
            // https://docs.seqan.de/seqan3/3-master-user/group__search__views.html#ga6e598d6a021868f704d39df73252974f
            it = hash;
        }
    };

    seqan::hibf::config config{.input_fn = my_input};
}
//![main]
