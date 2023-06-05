// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <type_traits>

#include <hibf/platform.hpp>

#include <cereal/details/traits.hpp>

namespace hibf
{

template <typename t>
concept cereal_output_archive = std::is_base_of_v<cereal::detail::OutputArchiveBase, t>;

template <typename t>
concept cereal_input_archive = std::is_base_of_v<cereal::detail::InputArchiveBase, t>;

template <typename t>
concept cereal_archive = cereal_output_archive<t> || cereal_input_archive<t>;

} // namespace hibf
