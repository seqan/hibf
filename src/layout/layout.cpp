// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>   // for find
#include <cassert>     // for assert
#include <charconv>    // for from_chars, from_chars_result
#include <cstddef>     // for size_t
#include <iostream>    // for operator<<, char_traits, basic_ostream, basic_istream, istream, ostream
#include <string>      // for basic_string, getline, string
#include <string_view> // for operator<<, string_view, operator==, basic_string_view
#include <vector>      // for vector

#include <hibf/layout/layout.hpp>   // for layout, operator<<
#include <hibf/layout/prefixes.hpp> // for layout_lower_level, layout_column_names, layout_fullest_technical_bin_idx

namespace seqan::hibf::layout
{

seqan::hibf::layout::layout::user_bin parse_layout_line(std::string const & current_line)
{
    seqan::hibf::layout::layout::user_bin result{};

    size_t tmp{}; // integer buffer when reading numbers

    // initialize parsing
    std::string_view const buffer{current_line};
    auto const buffer_end{buffer.end()};
    auto field_end = buffer.begin();
    assert(field_end != buffer_end);

    // read user bin index
    field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
    result.idx = tmp;
    assert(field_end != buffer_end && *field_end == '\t');

    do // read bin_indices
    {
        ++field_end; // skip tab or ;
        assert(field_end != buffer_end && *field_end != '\t');
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.previous_TB_indices.push_back(tmp);
    }
    while (field_end != buffer_end && *field_end != '\t');

    result.storage_TB_id = result.previous_TB_indices.back();
    result.previous_TB_indices.pop_back();

    do // read number of technical bins
    {
        ++field_end; // skip tab or ;
        field_end = std::from_chars(field_end, buffer_end, tmp).ptr;
        result.number_of_technical_bins = tmp; // only the last number really counts
    }
    while (field_end != buffer_end && *field_end != '\t');

    return result;
}

void seqan::hibf::layout::layout::read_from(std::istream & stream)
{
    // parse header
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

    std::getline(stream, line); // get first line that is always the max bin index of the top level bin
    assert(line.starts_with(prefix::layout_first_header_line));

    // parse High Level max bin index
    constexpr size_t fullest_tbx_prefix_size = prefix::layout_fullest_technical_bin_idx.size();
    assert(line.substr(prefix::layout_top_level.size() + 2, fullest_tbx_prefix_size)
           == prefix::layout_fullest_technical_bin_idx);
    std::string_view const hibf_max_bin_str{line.begin() + prefix::layout_top_level.size() + 2
                                                + fullest_tbx_prefix_size,
                                            line.end()};
    top_level_max_bin_id = parse_first_bin(hibf_max_bin_str);

    // read and parse header records, in order to sort them before adding them to the graph
    while (std::getline(stream, line) && line != prefix::layout_column_names)
    {
        assert(line.substr(1, prefix::layout_lower_level.size()) == prefix::layout_lower_level);

        // parse header line
        std::string_view const indices_str{
            line.begin() + 1 /*#*/ + prefix::layout_lower_level.size() + 1 /*_*/,
            std::find(line.begin() + prefix::layout_lower_level.size() + 2, line.end(), ' ')};

        assert(line.substr(prefix::layout_lower_level.size() + indices_str.size() + 3, fullest_tbx_prefix_size)
               == prefix::layout_fullest_technical_bin_idx);
        std::string_view const max_id_str{line.begin() + prefix::layout_lower_level.size() + indices_str.size()
                                              + fullest_tbx_prefix_size + 3,
                                          line.end()};

        max_bins.emplace_back(parse_bin_indices(indices_str), parse_first_bin(max_id_str));
    }

    assert(line == prefix::layout_column_names);

    // parse the rest of the file until either
    // 1) the end of the file is reached
    // 2) Another header line starts, which indicates a partitioned layout
    while (stream.good() && static_cast<char>(stream.peek()) != prefix::layout_header[0] &&  std::getline(stream, line))
        user_bins.emplace_back(parse_layout_line(line));
}

void seqan::hibf::layout::layout::write_to(std::ostream & stream) const
{
    // write layout header with max bin ids
    stream << prefix::layout_first_header_line << " " << prefix::layout_fullest_technical_bin_idx
           << top_level_max_bin_id << '\n';
    for (auto const & max_bin : max_bins)
        stream << max_bin << '\n';

    // write header line
    stream << prefix::layout_column_names << '\n';

    // write layout entries
    for (auto const & user_bin : user_bins)
        stream << user_bin << '\n';
}

void seqan::hibf::layout::layout::clear()
{
    top_level_max_bin_id = 0;
    max_bins.clear();
    user_bins.clear();
}

} // namespace seqan::hibf::layout
