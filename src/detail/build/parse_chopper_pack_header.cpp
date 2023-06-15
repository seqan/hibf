// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <algorithm>   // for find
#include <cassert>     // for assert
#include <charconv>    // for from_chars, from_chars_result
#include <cstddef>     // for size_t
#include <fstream>     // for basic_istream, istream
#include <string>      // for basic_string, char_traits, getline, operator==
#include <string_view> // for string_view, operator==, basic_string_view
#include <vector>      // for vector

#include <hibf/detail/build/parse_chopper_pack_header.hpp> // for parse_chopper_pack_header
#include <hibf/detail/layout/layout.hpp>                   // for layout
#include <hibf/detail/prefixes.hpp>                        // for merged_bin, header, header_config

namespace hibf
{

void parse_chopper_pack_header(std::istream & chopper_pack_file, hibf::layout::layout & hibf_layout)
{
    auto parse_bin_indices = [](std::string_view const & buffer)
    {
        std::vector<size_t> result;

        auto buffer_start = &buffer[0];
        auto const buffer_end = buffer_start + buffer.size();

        size_t tmp{};

        while (buffer_start < buffer_end)
        {
            buffer_start = std::from_chars(buffer_start, buffer_end, tmp).ptr;
            ++buffer_start; // skip ;
            result.push_back(tmp);
        }

        return result;
    };

    auto parse_first_bin = [](std::string_view const & buffer)
    {
        size_t tmp{};
        std::from_chars(&buffer[0], &buffer[0] + buffer.size(), tmp);
        return tmp;
    };

    std::string line;

    while (std::getline(chopper_pack_file, line) && line.size() >= 2
           && std::string_view{line}.substr(0, 1) == prefix::header
           && std::string_view{line}.substr(1, 1) == prefix::header_config)
        ; // skip config in header

    assert(line[0] == '#'); // we are reading header lines
    assert(line.substr(1, prefix::high_level.size()) == prefix::high_level);

    // parse High Level max bin index
    assert(line.substr(prefix::high_level.size() + 2, 11) == "max_bin_id:");
    std::string_view const hibf_max_bin_str{line.begin() + 27, line.end()};
    hibf_layout.top_level_max_bin_id = parse_first_bin(hibf_max_bin_str);

    // first read and parse header records, in order to sort them before adding them to the graph
    while (std::getline(chopper_pack_file, line) && line.substr(0, 6) != "#FILES")
    {
        assert(line.substr(1, prefix::merged_bin.size()) == prefix::merged_bin);

        // parse header line
        std::string_view const indices_str{line.begin() + 1 /*#*/ + prefix::merged_bin.size() + 1 /*_*/,
                                           std::find(line.begin() + prefix::merged_bin.size() + 2, line.end(), ' ')};

        assert(line.substr(prefix::merged_bin.size() + indices_str.size() + 3, 11) == "max_bin_id:");
        std::string_view const max_id_str{line.begin() + prefix::merged_bin.size() + indices_str.size() + 14,
                                          line.end()};

        hibf_layout.max_bins.emplace_back(parse_bin_indices(indices_str), parse_first_bin(max_id_str));
    }
}

} // namespace hibf
