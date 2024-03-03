// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cstdint>
#include <filesystem>

#include <sharg/all.hpp>

#include <seqan3/search/views/kmer_hash.hpp>

#include "index.hpp"
#include "sequence_io.hpp"

struct search_config
{
    std::filesystem::path index{};
    std::filesystem::path query{};
    std::filesystem::path output{};
    double threshold{0.7};
};

void search_index(search_config const & config)
{
    myindex index{};
    index.load(config.index);

    auto agent = index.hibf.membership_agent();

    id_seq_reader fin{config.query};
    std::vector<uint64_t> hashes{};
    std::ofstream fout{config.output};

    for (size_t i = 0; i < index.input_files.size(); ++i)
        fout << '#' << i << '\t' << index.input_files[i].c_str() << '\n';
    fout << "#QUERY_ID\tUSER_BINS\n";

    for (auto && [id, seq] : fin)
    {
        auto kmer_hash_view = seq | seqan3::views::kmer_hash(seqan3::ungapped{index.kmer});
        hashes.assign(kmer_hash_view.begin(), kmer_hash_view.end());
        uint16_t const threshold = static_cast<uint16_t>(config.threshold * hashes.size());

        auto & result = agent.membership_for(hashes, threshold);
        agent.sort_results(); // Optional

        fout << id << '\t';
        for (auto && r : result)
            fout << r << ' ';
        fout << '\n';
    }
}

void search(sharg::parser & parser)
{
    search_config config{};
    parser.add_option(config.index,
                      sharg::config{.short_id = 'i',
                                    .long_id = "index",
                                    .description = "Index",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(config.query,
                      sharg::config{.short_id = 'q',
                                    .long_id = "query",
                                    .description = "Query",
                                    .required = true,
                                    .validator = sharg::input_file_validator{{"fq", "fastq"}}});

    parser.add_option(config.output,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "Output",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});

    parser.add_option(config.threshold,
                      sharg::config{.short_id = '\0',
                                    .long_id = "threshold",
                                    .description = "Threshold.",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "Parsing error. " << ext.what() << '\n';
        std::exit(EXIT_FAILURE);
    }

    search_index(config);
}
