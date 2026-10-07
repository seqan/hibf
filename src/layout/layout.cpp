// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>   // for find_if, sort, adjacent_find, equal_range, mismatch, max
#include <cassert>     // for assert
#include <charconv>    // for from_chars, from_chars_result
#include <cmath>       // for floor
#include <compare>     // for operator<=>
#include <cstddef>     // for size_t
#include <format>      // for format, format_to
#include <iostream>    // for operator<<, char_traits, basic_ostream, basic_istream, istream, ostream, cerr
#include <iterator>    // for back_inserter
#include <limits>      // for numeric_limits
#include <map>         // for map
#include <memory>      // for addressof
#include <optional>    // for optional, nullopt
#include <ranges>      // for iota, to, transform
#include <set>         // for set
#include <stdexcept>   // for invalid_argument
#include <string>      // for basic_string, getline, string
#include <string_view> // for operator<<, string_view, operator==, basic_string_view
#include <utility>     // for move
#include <vector>      // for vector

#include <hibf/config.hpp>                                // for config
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for bin_kind
#include <hibf/layout/layout.hpp>                         // for layout, operator<<
#include <hibf/layout/prefixes.hpp>          // for layout_lower_level, layout_column_names, layout_fullest_techni...
#include <hibf/misc/add_empty_bins.hpp>      // for add_empty_bins
#include <hibf/misc/next_multiple_of_64.hpp> // for next_multiple_of_64
#include <hibf/misc/subtract_empty_bins.hpp> // for subtract_empty_bins

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

        auto buffer_start = buffer.data();
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
        std::from_chars(buffer.data(), buffer.data() + buffer.size(), tmp);
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
    while (stream.good() && static_cast<char>(stream.peek()) != prefix::layout_header[0] && std::getline(stream, line))
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

namespace
{

// A range of technical bins within one IBF: a merged bin, or the technical bins of a single or split user bin.
// Like hierarchical_interleaved_bloom_filter::ibf_bin_to_user_bin_id, merged bins have user bin bin_kind::merged.
struct occupied_bins
{
    size_t first{};
    size_t last{};                     // inclusive
    size_t user_bin{bin_kind::merged}; // a user bin index or bin_kind::merged (std::numeric_limits<uint64_t>::max())

    bool is_merged() const
    {
        return user_bin == bin_kind::merged;
    }

    friend auto operator<=>(occupied_bins const &, occupied_bins const &) = default;
};

// The merged bin or user bin starting at `technical_bin` in `bins` (sorted, non-overlapping), or nullptr.
occupied_bins const * find_first_occupied_bin(std::vector<occupied_bins> const & bins, size_t const technical_bin)
{
    auto const range = std::ranges::equal_range(bins, technical_bin, {}, &occupied_bins::first);
    return range.empty() ? nullptr : std::addressof(range.front());
}

// "the Root-IBF" for the top-level IBF ("The Root-IBF" at the start of a sentence), or "IBF 2;3" for the IBF below
// merged bin 3 of the IBF below merged bin 2 of the top-level IBF. "IBF 0" is the IBF below merged bin 0 of the
// top-level IBF.
std::string ibf_name(std::vector<size_t> const & path, bool const sentence_start = false)
{
    std::string result{path.empty() ? (sentence_start ? "The Root-IBF" : "the Root-IBF") : "IBF "};
    for (size_t const technical_bin : path)
        std::format_to(std::back_inserter(result), "{};", technical_bin);
    if (!path.empty())
        result.pop_back();
    return result;
}

} // namespace

bool seqan::hibf::layout::layout::validate(config const & config, diagnostic_handler const & handler) const
{
    using code = diagnostic::code;
    using severity = diagnostic::severity;

    auto report = [&handler](diagnostic const & finding)
    {
        if (handler)
            handler(finding);
    };

    auto error = [&report](diagnostic finding)
    {
        finding.level = severity::error;
        report(finding);
        return false;
    };

    // The largest number of technical bins of an IBF that next_multiple_of_64 can round up without overflow.
    constexpr size_t max_ibf_technical_bins{std::numeric_limits<size_t>::max() - 63u};

    if (user_bins.empty())
        return error({.what = code::empty_layout, .message = "The layout contains no user bins."});

    // Checks that only concern a single user bin.
    for (auto const & ub : user_bins)
    {
        if (ub.number_of_technical_bins == 0u)
            return error({.what = code::no_technical_bins,
                          .ibf = ub.previous_TB_indices,
                          .user_bin = ub.idx,
                          .technical_bin = ub.storage_TB_id,
                          .message = std::format("User bin {} occupies zero technical bins.", ub.idx)});

        // The graph computes the number of technical bins of an IBF as `storage_TB_id + number_of_technical_bins` and
        // `merged bin + 1`. The build rounds it up to a multiple of 64.
        if (ub.storage_TB_id >= max_ibf_technical_bins
            || ub.number_of_technical_bins > max_ibf_technical_bins - ub.storage_TB_id)
            return error({.what = code::technical_bin_overflow,
                          .ibf = ub.previous_TB_indices,
                          .user_bin = ub.idx,
                          .technical_bin = ub.storage_TB_id,
                          .message = std::format("In {}, the technical bins of user bin {} exceed the maximum number "
                                                 "of technical bins of an IBF.",
                                                 ibf_name(ub.previous_TB_indices),
                                                 ub.idx)});

        if (auto it = std::ranges::find_if(ub.previous_TB_indices,
                                           [](size_t const merged_bin)
                                           {
                                               return merged_bin >= max_ibf_technical_bins;
                                           });
            it != ub.previous_TB_indices.end())
        {
            std::vector<size_t> path(ub.previous_TB_indices.begin(), it);
            std::string const name = ibf_name(path);
            return error({.what = code::technical_bin_overflow,
                          .ibf = std::move(path),
                          .user_bin = ub.idx,
                          .technical_bin = *it,
                          .message = std::format("In {}, the merged bin on the path of user bin {} exceeds the "
                                                 "maximum number of technical bins of an IBF.",
                                                 name,
                                                 ub.idx)});
        }

        if (ub.idx >= config.number_of_user_bins)
            return error({.what = code::user_bin_out_of_range,
                          .ibf = ub.previous_TB_indices,
                          .user_bin = ub.idx,
                          .message = std::format("User bin index {} is not below the number of user bins ({}).",
                                                 ub.idx,
                                                 config.number_of_user_bins)});
    }

    // Checks that concern the set of user bin indices.
    {
        auto indices = user_bins | std::views::transform(&user_bin::idx) | std::ranges::to<std::vector>();
        std::ranges::sort(indices);

        if (auto it = std::ranges::adjacent_find(indices); it != indices.end())
            return error({.what = code::duplicate_user_bin,
                          .user_bin = *it,
                          .message = std::format("User bin {} occurs more than once.", *it)});

        if (indices.size() != config.number_of_user_bins)
        {
            // The indices are unique and below number_of_user_bins. The first index that differs from its position
            // is missing.
            size_t const missing = std::ranges::mismatch(indices, std::views::iota(size_t{})).in1 - indices.begin();
            return error({.what = code::missing_user_bin,
                          .user_bin = missing,
                          .message = std::format("User bin {} is missing from the layout.", missing)});
        }
    }

    // An IBF is identified by the merged bins on the **path** from the top-level IBF.
    // The same path is stored for each user bin in layout::user_bin::previous_TB_indices.
    std::map<std::vector<size_t>, std::vector<occupied_bins>> ibfs{};

    // Fill map `ibfs` by adding `occupied_bins` to each IBF identified by its `path`.
    for (auto const & ub : user_bins)
    {
        std::vector<size_t> path{};
        auto ibf = ibfs.try_emplace(path).first;
        for (size_t const technical_bin : ub.previous_TB_indices)
        {
            path.push_back(technical_bin);
            // Every user bin below a merged bin passes through it. Only the first one creates the IBF below the merged
            // bin, and only then is the merged bin added. Hence, each merged bin occurs once.
            auto const [below, created] = ibfs.try_emplace(path);
            if (created)
                ibf->second.push_back({technical_bin, technical_bin, bin_kind::merged});
            ibf = below;
        }

        size_t const last = ub.storage_TB_id + ub.number_of_technical_bins - 1u;
        ibf->second.push_back({ub.storage_TB_id, last, ub.idx});
    }

    // Within each IBF, the occupied ranges must not overlap.
    // If all adjacent ranges (sorted by first technical bin) are disjoint, all ranges are disjoint.
    for (auto & [path, bins] : ibfs)
    {
        std::ranges::sort(bins);
        for (size_t i = 1u; i < bins.size(); ++i)
        {
            occupied_bins const & previous = bins[i - 1u];
            occupied_bins const & current = bins[i];

            if (current.first <= previous.last)
            {
                // At most one of them is a merged bin: Two distinct merged bins cannot overlap.
                auto describe = [](occupied_bins const & bin)
                {
                    return bin.is_merged() ? std::string{"a merged bin"} : std::format("user bin {}", bin.user_bin);
                };
                return error({.what = code::overlapping_technical_bins,
                              .ibf = path,
                              .user_bin = current.is_merged() ? previous.user_bin : current.user_bin,
                              .technical_bin = current.first,
                              .message = std::format("In {}, technical bin {} is used by {} and by {}.",
                                                     ibf_name(path),
                                                     current.first,
                                                     describe(previous),
                                                     describe(current))});
            }
        }
    }

    // The build sizes its FPR correction table by the number of technical bins of the Root-IBF and looks it up with the
    // number of technical bins of each max bin (build_index, construct_ibf). This only becomes relevant once we expect
    // layouts where tmax is not a strict upper bound for the number of technical bins of an IBF.
    size_t const root_technical_bins = ibfs.at({}).back().last + 1u;       // bins are sorted and disjoint
    size_t const max_bin_limit = next_multiple_of_64(root_technical_bins); // no overflow, see technical_bin_overflow

    auto invalid_max_bin = [&error](std::vector<size_t> const & path, size_t const id)
    {
        return error({.what = code::invalid_max_bin,
                      .ibf = path,
                      .technical_bin = id,
                      .message = std::format("The max bin (\"{}\") of {} is neither a merged bin nor the first "
                                             "technical bin of a user bin.",
                                             prefix::layout_fullest_technical_bin_idx,
                                             ibf_name(path))});
    };

    std::set<std::vector<size_t>> ibfs_with_max_bin{};
    for (auto const & [path, id] : max_bins)
    {
        if (path.empty())
            return error({.what = code::top_level_max_bin_entry,
                          .ibf = std::vector<size_t>{},
                          .technical_bin = id,
                          .message = std::format("The max bins contain an entry for the Root-IBF. The Root-IBF's max "
                                                 "bin is given by top_level_max_bin_id (\"{}\"), not by the max bins "
                                                 "(\"{}{}\").",
                                                 prefix::layout_first_header_line,
                                                 prefix::layout_header,
                                                 prefix::layout_lower_level)});

        auto it = ibfs.find(path);
        if (it == ibfs.end())
            return error({.what = code::max_bin_without_ibf,
                          .ibf = path,
                          .technical_bin = id,
                          .message = std::format("The max bins (\"{}{}\") contain an entry for {}, which does not "
                                                 "exist.",
                                                 prefix::layout_header,
                                                 prefix::layout_lower_level,
                                                 ibf_name(path))});

        if (!ibfs_with_max_bin.insert(path).second)
            return error({.what = code::duplicate_max_bin,
                          .ibf = path,
                          .technical_bin = id,
                          .message = std::format("The max bins (\"{}{}\") contain more than one entry for {}.",
                                                 prefix::layout_header,
                                                 prefix::layout_lower_level,
                                                 ibf_name(path))});

        occupied_bins const * const max_bin = find_first_occupied_bin(it->second, id);
        if (max_bin == nullptr)
            return invalid_max_bin(path, id);

        if (size_t const span = max_bin->last - max_bin->first + 1u; span > max_bin_limit)
            return error({.what = code::max_bin_exceeds_root,
                          .ibf = path,
                          .user_bin = max_bin->user_bin,
                          .technical_bin = id,
                          .message = std::format("The max bin (\"{}\") of {} spans {} technical bins, but at most {} "
                                                 "are supported (the Root-IBF's {} technical bins, rounded up to a "
                                                 "multiple of 64).",
                                                 prefix::layout_fullest_technical_bin_idx,
                                                 ibf_name(path),
                                                 span,
                                                 max_bin_limit,
                                                 root_technical_bins)});
    }

    for (auto const & [path, bins] : ibfs)
        if (!path.empty() && !ibfs_with_max_bin.contains(path))
            return error({.what = code::missing_max_bin,
                          .ibf = path,
                          .message = std::format("The max bins (\"{}{}\") contain no entry for {}.",
                                                 prefix::layout_header,
                                                 prefix::layout_lower_level,
                                                 ibf_name(path))});

    if (size_t const id = top_level_max_bin_id; find_first_occupied_bin(ibfs.at({}), id) == nullptr)
        return invalid_max_bin({}, id);

    // Warnings and notes do not affect the result.
    if (!handler)
        return true;

    // The lowest levels use up to next_multiple_of_64(#user bins) technical bins, which is only bounded by tmax if
    // tmax is a multiple of 64. config::validate_and_set_defaults() rounds tmax up to a multiple of 64.
    bool const check_tmax = config.tmax != 0u && config.tmax <= max_ibf_technical_bins;
    size_t const max_technical_bins =
        check_tmax ? next_multiple_of_64(config.tmax) : std::numeric_limits<size_t>::max();
    std::string const tmax_name = [&]()
    {
        std::string result{"tmax"};
        if (check_tmax && config.tmax != max_technical_bins)
            result = std::format("tmax {} (rounded up to a multiple of 64)", config.tmax);
        return result;
    }();
    // The build adds empty bins to each IBF (interleaved_bloom_filter's constructor). An IBF whose last used technical
    // bin is `last` has next_multiple_of_64(add_empty_bins(last + 1, empty_bin_fraction)) technical bins.
    // std::nullopt if this does not fit into size_t: add_empty_bins is checked in floating point first.
    double const empty_bin_fraction =
        (config.empty_bin_fraction > 0.0 && config.empty_bin_fraction < 1.0) ? config.empty_bin_fraction : 0.0;
    auto const built_technical_bins = [empty_bin_fraction](size_t const last) -> std::optional<size_t>
    {
        if (empty_bin_fraction == 0.0)
            return next_multiple_of_64(last + 1u); // no overflow, see technical_bin_overflow
        if (std::floor((last + 1u) / (1.0 - empty_bin_fraction)) >= 0x1p63)
            return std::nullopt;
        return next_multiple_of_64(add_empty_bins(last + 1u, empty_bin_fraction));
    };

    for (auto const & [path, bins] : ibfs)
    {
        std::string const name = ibf_name(path, true); // all messages below start with it
        size_t const last = bins.back().last;          // bins are sorted and disjoint
        std::optional<size_t> const built = built_technical_bins(last);

        // max_technical_bins is a multiple of 64: last >= max_technical_bins implies *built > max_technical_bins.
        if (check_tmax && (!built.has_value() || *built > max_technical_bins))
            report({.level = severity::warning,
                    .what = code::technical_bin_exceeds_tmax,
                    .ibf = path,
                    .technical_bin = last,
                    .message = (last >= max_technical_bins)
                                 ? std::format("{} uses technical bin {}, but {} allows only {} technical bins.",
                                               name,
                                               last,
                                               tmax_name,
                                               max_technical_bins)
                                 : std::format("{} uses technical bin {}. With the empty bins added by the build "
                                               "(empty_bin_fraction {}), it exceeds the {} technical bins {} allows.",
                                               name,
                                               last,
                                               empty_bin_fraction,
                                               max_technical_bins,
                                               tmax_name)});

        if (!path.empty() && bins.size() == 1u)
        {
            occupied_bins const & bin = bins.front();
            report({.level = severity::warning,
                    .what = code::single_bin_ibf,
                    .ibf = path,
                    .user_bin = bin.is_merged() ? std::nullopt : std::optional<size_t>{bin.user_bin},
                    .technical_bin = bin.first,
                    .message = bin.is_merged()
                                 ? std::format("{} only contains merged bin {}. Moving the IBF below it up would "
                                               "save a level.",
                                               name,
                                               bin.first)
                                 : std::format("{} only contains user bin {}. Storing it in the parent IBF would "
                                               "save a level.",
                                               name,
                                               bin.user_bin)});
        }

        // Empty technical bins below the last used technical bin. The used technical bins and the empty bins should
        // partition the technical bins of an IBF.
        size_t used{};
        std::optional<size_t> first_empty{};
        for (size_t i = 0u; i < bins.size(); ++i)
        {
            size_t const next = (i == 0u) ? 0u : bins[i - 1u].last + 1u;
            if (!first_empty.has_value() && bins[i].first > next)
                first_empty = next;
            used += bins[i].last - bins[i].first + 1u;
        }
        size_t const intermittent_empty = last + 1u - used;

        // The empty technical bins after the last used technical bin must be as many as empty_bin_fraction implies,
        // i.e., none if it is 0. Intermittent empty technical bins do not count.
        if (built.has_value())
        {
            size_t const trailing_empty = *built - (last + 1u);
            size_t const expected_empty = *built - subtract_empty_bins(*built, empty_bin_fraction);
            if (trailing_empty != expected_empty)
                report({.level = severity::warning,
                        .what = code::unexpected_empty_bins,
                        .ibf = path,
                        .technical_bin = last,
                        .message = std::format("{} ends with {} empty technical bin{}, but {} {} expected. There is a "
                                               "total of {} technical bins and the empty_bin_fraction is {}.",
                                               name,
                                               trailing_empty,
                                               (trailing_empty == 1u) ? "" : "s",
                                               expected_empty,
                                               (expected_empty == 1u) ? "is" : "are",
                                               *built,
                                               empty_bin_fraction)});
        }

        if (first_empty.has_value())
            report({.level = severity::note,
                    .what = code::empty_technical_bins,
                    .ibf = path,
                    .technical_bin = first_empty,
                    .message = (intermittent_empty == 1u)
                                 ? std::format("{} uses technical bins 0-{}, but technical bin {} is empty.",
                                               name,
                                               last,
                                               *first_empty)
                                 : std::format("{} uses technical bins 0-{}, but {} of them are empty; the first one "
                                               "is {}.",
                                               name,
                                               last,
                                               intermittent_empty,
                                               *first_empty)});
    }

    return true;
}

size_t seqan::hibf::layout::layout::number_of_levels() const
{
    size_t levels{};
    for (auto const & ub : user_bins)
        levels = std::max(levels, ub.previous_TB_indices.size() + 1u);
    return levels;
}

void seqan::hibf::layout::layout::throw_on_error(diagnostic const & finding)
{
    if (finding.level == diagnostic::severity::error)
        throw std::invalid_argument{std::format("{}", finding)};
    if (finding.level == diagnostic::severity::warning)
        std::cerr << finding << '\n';
#ifndef NDEBUG
    if (finding.level == diagnostic::severity::note)
        std::cerr << finding << '\n';
#endif
}

} // namespace seqan::hibf::layout
