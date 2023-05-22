// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <lemon/list_graph.h> /// Must be first include.

#include <fstream>
#include <string>

#include <hibf/detail/build/hibf/parse_chopper_pack_header.hpp>
#include <hibf/detail/build/hibf/parse_chopper_pack_line.hpp>

namespace hibf
{

inline layout::layout read_chopper_pack_file(std::vector<std::vector<std::string>> & filenames,
                                             std::string const & chopper_pack_filename)
{
    layout::layout hibf_layout{};

    std::ifstream chopper_pack_file{chopper_pack_filename};

    if (!chopper_pack_file.good() || !chopper_pack_file.is_open())
        throw std::logic_error{"Could not open file " + chopper_pack_filename + " for reading"}; // GCOVR_EXCL_LINE

    // parse header
    // -------------------------------------------------------------------------
    parse_chopper_pack_header(chopper_pack_file, hibf_layout);

    std::string current_line;
    while (std::getline(chopper_pack_file, current_line))
        hibf_layout.user_bins.emplace_back(parse_chopper_pack_line(current_line, filenames));

    return hibf_layout;
}

} // namespace hibf
