// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <hibf/platform.hpp>

#include <cereal/details/helpers.hpp> // for InputArchiveBase, OutputArchiveBase
#include <cereal/details/traits.hpp>  // for InputArchiveBase, OutputArchiveBase

namespace seqan::hibf
{

template <typename t>
concept cereal_output_archive = std::is_base_of_v<cereal::detail::OutputArchiveBase, t>;

template <typename t>
concept cereal_input_archive = std::is_base_of_v<cereal::detail::InputArchiveBase, t>;

template <typename t>
concept cereal_archive = cereal_output_archive<t> || cereal_input_archive<t>;

template <typename t>
concept cereal_text_archive = std::is_base_of_v<cereal::traits::TextArchive, t>;

} // namespace seqan::hibf
