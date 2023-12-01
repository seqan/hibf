// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <concepts> // for derived_from
#include <cstddef>  // for size_t
#include <iosfwd>   // for ostream, istream
#include <vector>   // for operator==, vector

#include <hibf/layout/prefixes.hpp> // for layout_fullest_technical_bin_idx, layout_header, layout_lower_level
#include <hibf/platform.hpp>

namespace seqan::hibf::layout
{

/*!\brief The layout.
 * \ingroup hibf_layout
 */
struct layout
{
    struct max_bin
    {
        std::vector<size_t> previous_TB_indices{}; // identifies the IBF based on upper levels
        size_t id{};                               // the technical bin id that has the maximum kmer content

        friend auto operator<=>(max_bin const &, max_bin const &) = default;

        // needs a template (instead of using std::ostream directly) to be able to only include <iosfwd>
        template <typename stream_type>
            requires std::derived_from<stream_type, std::ostream>
        friend stream_type & operator<<(stream_type & stream, max_bin const & object)
        {
            stream << prefix::layout_header << prefix::layout_lower_level << '_';
            auto it = object.previous_TB_indices.begin();
            auto end = object.previous_TB_indices.end();
            // If not empty, we join with ';'
            if (it != end)
            {
                stream << *it;
                while (++it != end)
                    stream << ';' << *it;
            }
            stream << " " << prefix::layout_fullest_technical_bin_idx << object.id;

            return stream;
        }
    };

    struct user_bin
    {
        std::vector<size_t> previous_TB_indices{}; // previous technical bin indices which refer to merged bin indices.
        size_t storage_TB_id{};                    // id of the technical bin that the user bin is actuallly stored in
        size_t number_of_technical_bins{};         // 1 == single bin, >1 == split_bin
        size_t idx{};                              // The index of the user bin corresponding to the order in data

        friend auto operator<=>(user_bin const &, user_bin const &) = default;

        // needs a template (instead of using std::ostream directly) to be able to only include <iosfwd>
        template <typename stream_type>
            requires std::derived_from<stream_type, std::ostream>
        friend stream_type & operator<<(stream_type & stream, user_bin const & object)
        {
            stream << object.idx << '\t';
            for (auto bin : object.previous_TB_indices)
                stream << bin << ';';
            stream << object.storage_TB_id << '\t';
            for ([[maybe_unused]] auto && elem : object.previous_TB_indices) // number of bins per merged level is 1
                stream << "1;";
            stream << object.number_of_technical_bins;

            return stream;
        }
    };

    void read_from(std::istream & stream);
    void write_to(std::ostream & stream) const;

    void clear();

    size_t top_level_max_bin_id{};
    std::vector<max_bin> max_bins{};
    std::vector<user_bin> user_bins{};

    bool operator==(layout const &) const = default;
};

} // namespace seqan::hibf::layout
