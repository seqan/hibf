// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <sharg/parser.hpp>

#include <cereal/archives/binary.hpp>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

#include "fmt.hpp"
#include "inspect.hpp"

struct config
{
    std::filesystem::path input{};
    bool is_hibf{};
};

config parse_arguments(std::vector<std::string> command_line)
{
    sharg::parser parser{"inspect", std::move(command_line), sharg::update_notifications::off};
    config cfg{};

    parser.add_option(cfg.input,
                      sharg::config{.short_id = '\0',
                                    .long_id = "input",
                                    .description = "The index to inspect.",
                                    .validator = sharg::input_file_validator{}});

    parser.add_flag(cfg.is_hibf,
                    sharg::config{.short_id = '\0', .long_id = "hibf", .description = "The index is an HIBF"});

    parser.info.author = "Enrico Seiler";
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.short_description = "Inspect an IBF/HIBF.";
    parser.parse();

    return cfg;
}

int main(int argc, char ** argv)
{
    try
    {
        config const config = parse_arguments({argv, argv + argc});

        if (!config.is_hibf)
        {
            seqan::hibf::interleaved_bloom_filter ibf{};
            {
                std::ifstream os{config.input, std::ios::binary};
                cereal::BinaryInputArchive iarchive{os};
                iarchive(ibf);
            }

            seqan::hibf::util::inspect(ibf);
        }
        else
        {
            seqan::hibf::hierarchical_interleaved_bloom_filter hibf{};
            {
                std::ifstream os{config.input, std::ios::binary};
                cereal::BinaryInputArchive iarchive{os};
                iarchive(hibf);
            }

            seqan::hibf::util::inspect(hibf);
        }
    }
    catch (std::exception const & ext)
    {
        fmt::print(stderr, seqan::hibf::util::styles::color(fmt::color::red, stderr), "[Error]");
        fmt::print(stderr, " {}\n", ext.what());
        std::exit(-1);
    }

    return 0;
}
