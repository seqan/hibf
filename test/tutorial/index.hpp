// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <hibf/cereal/path.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

class myindex
{
public:
    uint8_t kmer{};
    std::vector<std::filesystem::path> input_files{};
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{};

    myindex() = default;
    myindex & operator=(myindex const &) = default;
    myindex(myindex const &) = default;
    myindex(myindex &&) = default;
    myindex & operator=(myindex &&) = default;
    ~myindex() = default;

    explicit myindex(uint8_t const kmer,
                     std::vector<std::filesystem::path> input_files,
                     seqan::hibf::hierarchical_interleaved_bloom_filter hibf) :
        kmer{kmer},
        input_files{std::move(input_files)},
        hibf{std::move(hibf)}
    {}

    void store(std::filesystem::path const & path) const
    {
        std::ofstream fout{path};
        cereal::BinaryOutputArchive oarchive{fout};
        oarchive(*this);
    }

    void load(std::filesystem::path const & path)
    {
        std::ifstream fin{path};
        cereal::BinaryInputArchive iarchive{fin};
        iarchive(*this);
    }

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(kmer);
        archive(input_files);
        archive(hibf);
    }
};
