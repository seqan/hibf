// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <numeric>

#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

namespace seqan::hibf
{

struct inspector
{
    using ibf_t = seqan::hibf::interleaved_bloom_filter;
    using hibf_t = seqan::hibf::hierarchical_interleaved_bloom_filter;

    static size_t const hash_function_count(ibf_t const & ibf)
    {
        return ibf.hash_function_count();
    }

    static size_t const bin_count(ibf_t const & ibf)
    {
        return ibf.bin_count();
    }

    static size_t const bin_size(ibf_t const & ibf)
    {
        return ibf.bin_size();
    }

    static size_t const bit_size(ibf_t const & ibf)
    {
        return ibf.bit_size();
    }

    static size_t const technical_bins(ibf_t const & ibf)
    {
        return ibf.technical_bins;
    }

    static size_t const empty_bins(ibf_t const & ibf)
    {
        return ibf.technical_bins - ibf.bin_count();
    }

    static size_t const hash_shift(ibf_t const & ibf)
    {
        return ibf.hash_shift;
    }

    static size_t const bin_words(ibf_t const & ibf)
    {
        return ibf.bin_words;
    }

    static constexpr std::array<size_t, 5> const & hash_seeds(ibf_t const & ibf)
    {
        return ibf.hash_seeds;
    }

    static std::vector<size_t> const & occupancy(ibf_t const & ibf)
    {
        return ibf.occupancy;
    }

    static size_t const number_of_user_bins(hibf_t const & hibf)
    {
        return hibf.number_of_user_bins;
    }

    static size_t const number_of_ibfs(hibf_t const & hibf)
    {
        return hibf.ibf_vector.size();
    }

    static size_t const total_size(hibf_t const & hibf)
    {
        auto vec_size = [](auto const & vec)
        {
            using value_t = std::ranges::range_value_t<decltype(vec)>;
            return vec.size() * sizeof(value_t);
        };

        size_t const ibf_vector_size = std::transform_reduce(hibf.ibf_vector.cbegin(), //
                                                             hibf.ibf_vector.cend(),
                                                             size_t{},
                                                             std::plus{},
                                                             bit_size);

        size_t const next_ibf_id_size = std::transform_reduce(hibf.next_ibf_id.cbegin(), //
                                                              hibf.next_ibf_id.cend(),
                                                              size_t{},
                                                              std::plus{},
                                                              vec_size);

        size_t const prev_ibf_id_size = vec_size(hibf.prev_ibf_id);

        size_t const ibf_bin_to_user_bin_id_size = std::transform_reduce(hibf.ibf_bin_to_user_bin_id.cbegin(),
                                                                         hibf.ibf_bin_to_user_bin_id.cend(),
                                                                         size_t{},
                                                                         std::plus{},
                                                                         vec_size);

        return ibf_vector_size + next_ibf_id_size + prev_ibf_id_size + ibf_bin_to_user_bin_id_size;
    }

    static size_t const total_empty(hibf_t const & hibf)
    {
        return std::transform_reduce(hibf.ibf_vector.cbegin(),
                                     hibf.ibf_vector.cend(),
                                     size_t{},
                                     std::plus{},
                                     empty_bins);
    }

    static std::vector<ibf_t> const & ibf_vector(hibf_t const & hibf)
    {
        return hibf.ibf_vector;
    }

    static std::vector<std::vector<uint64_t>> const & next_ibf_id(hibf_t const & hibf)
    {
        return hibf.next_ibf_id;
    }

    static std::vector<hibf_t::previous_ibf_id_pair> const & prev_ibf_id(hibf_t const & hibf)
    {
        return hibf.prev_ibf_id;
    }

    static std::vector<std::vector<uint64_t>> const & ibf_bin_to_user_bin_id(hibf_t const & hibf)
    {
        return hibf.ibf_bin_to_user_bin_id;
    }

    static size_t const user_bins(hibf_t const & hibf, size_t const idx)
    {
        auto & vec = hibf.ibf_bin_to_user_bin_id[idx];
        size_t result{};
        for (size_t i = 0u; i < vec.size(); ++i)
        {
            bool const is_special = vec[i] == seqan::hibf::bin_kind::merged || vec[i] == seqan::hibf::bin_kind::deleted;
            bool const is_split = (i != vec.size() - 1u) && vec[i] == vec[i + 1u];
            result += !is_special && !is_split;
        }
        return result;
    }
};

} // namespace seqan::hibf
