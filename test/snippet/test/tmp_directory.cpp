#include <gtest/gtest.h> // for Test, AssertionResult, TestInfo, EXPECT_TRUE, Message, TEST, TestPar...

#include <filesystem> // for remove
#include <fstream>    // for char_traits, basic_ofstream, basic_ostream, operator<<, ofstream
#include <memory>     // for allocator

#include <hibf/test/sandboxed_path.hpp> // for operator/, sandboxed_path
#include <hibf/test/tmp_directory.hpp>  // for tmp_directory

TEST(snippet_tmp_directory, tmp_directory_)
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

    // check that everything was cleanup properly
    EXPECT_TRUE(tmp.empty());
}
