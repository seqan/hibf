// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <ranges>
#include <seqan3/std/charconv>

#include <seqan3/search/dream_index/detail/build/hibf/chopper_pack_record.hpp>

namespace seqan3::hibf
{

chopper_pack_record parse_chopper_pack_line(std::string const & current_line)
{
    chopper_pack_record result{};

    // initialize parsing
    std::string_view const buffer{current_line};
    auto const buffer_end{buffer.end()};
    auto field_end = buffer.begin();
    while (field_end != buffer_end && *field_end != '\t')
        ++field_end;

    // parse filenames
    std::string_view const filenames{buffer.begin(), field_end};
    for (auto const && filename : filenames | std::views::split(';'))
    {
        auto const common_view = filename | std::views::common;
        result.filenames.emplace_back(common_view.begin(), common_view.end());
    }

    size_t tmp{};

    // parse filename for user bin index
    assert(result.filenames.size() == 1u);
    assert(result.filenames[0][0] == 'U');
    assert(result.filenames[0][1] == 'B');
    assert(result.filenames[0][2] == '_');
    std::from_chars(result.filenames[0].data() + 3 /* UB_ */,
                    result.filenames[0].data() + result.filenames[0].size(),
                    tmp);
    result.user_bin_index = tmp;

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.bin_indices.push_back(tmp);
    }
    while (field_end != buffer_end && *field_end != '\t');

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.number_of_bins.push_back(tmp);
    }
    while (field_end != buffer_end && *field_end != '\t');

    return result;
}

} // namespace seqan3::hibf
