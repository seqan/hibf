// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cstdint>
#include <filesystem>

#include <sharg/all.hpp>

#include <seqan3/search/views/kmer_hash.hpp>

#include "index.hpp"
#include "sequence_io.hpp"
#include "validator.hpp"

struct build_config
{
    std::filesystem::path input{};
    std::vector<std::filesystem::path> input_files{};
    std::filesystem::path output{};
    uint8_t kmer{};
    uint8_t threads{1u};
};

void build_hibf(build_config & config)
{
    auto input_lambda = [&config](size_t const user_bin_index, seqan::hibf::insert_iterator it)
    {
        seq_reader fin{config.input_files[user_bin_index]};
        for (auto && [seq] : fin)
            for (auto && hash : seq | seqan3::views::kmer_hash(seqan3::ungapped{config.kmer}))
                it = hash;
    };

    seqan::hibf::config hibf_config{.input_fn = input_lambda,
                                    .number_of_user_bins = config.input_files.size(),
                                    .number_of_hash_functions = 2u,
                                    .maximum_fpr = 0.05,
                                    .threads = config.threads};

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{hibf_config};
    myindex index{config.kmer, std::move(config.input_files), std::move(hibf)};
    index.store(config.output);
}

void read_input_files(build_config & config)
{
    std::ifstream file{config.input};
    std::string line{};

    while (std::getline(file, line))
        config.input_files.emplace_back(line);

    if (config.input_files.empty())
        throw sharg::validation_error{"No input files found in " + config.input.string() + "."};

    std::ranges::for_each(config.input_files, sharg::input_file_validator{{"fa", "fasta"}});
}

void build(sharg::parser & parser)
{
    build_config config{};
    parser.add_option(config.input,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input",
                                    .description = "Input",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(config.output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "Output",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});

    parser.add_option(config.kmer,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer",
                                    .description = "Kmer",
                                    .required = true,
                                    .validator = sharg::arithmetic_range_validator{1, 32}});

    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Threads.",
                                    .validator = positive_integer_validator<decltype(config.threads)>{}});

    try
    {
        parser.parse();
        read_input_files(config);
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "Parsing error. " << ext.what() << '\n';
        std::exit(EXIT_FAILURE);
    }

    build_hibf(config);
}
