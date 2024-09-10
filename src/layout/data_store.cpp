// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <hibf/layout/data_store.hpp> // for data_store

namespace seqan::hibf::layout
{

void data_store::validate() const
{
    if (hibf_layout == nullptr)
        throw std::invalid_argument{"[HIBF ERROR] data_store::hibf_layout must not be nullptr."};

    if (kmer_counts == nullptr)
        throw std::invalid_argument{"[HIBF ERROR] data_store::kmer_counts must not be nullptr."};

    if (sketches != nullptr && kmer_counts->size() != sketches->size())
        throw std::invalid_argument{
            "[HIBF ERROR] data_store::kmer_counts and data_store::sketches must have the same size."};

    if (positions.size() > kmer_counts->size())
        throw std::invalid_argument{
            "[HIBF ERROR] data_store::kmer_counts.size() must not be smaller than data_store::positions.size()."};

    if (fpr_correction.empty())
        throw std::invalid_argument{"[HIBF ERROR] data_store::fpr_correction must not be empty."};

    if (relaxed_fpr_correction <= 0.0 || relaxed_fpr_correction > 1.0)
        throw std::invalid_argument{"[HIBF ERROR] data_store::relaxed_fpr_correction must be in (0.0,1.0]."};
}

} // namespace seqan::hibf::layout
