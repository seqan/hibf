#include <fstream>
#include <string>
#include <vector>

#include <hibf/test/temporary_snippet_file.hpp>

seqan::hibf::test::temporary_snippet_file file1{"file1.fa", "ACGT"};
seqan::hibf::test::temporary_snippet_file file2{"file2.fa", "ACGT"};

//![main]
#include <hibf/config.hpp>

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
