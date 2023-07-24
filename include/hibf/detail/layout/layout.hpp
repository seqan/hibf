#pragma once

#include <algorithm> // for copy
#include <concepts>  // for derived_from
#include <cstddef>   // for size_t
#include <iosfwd>    // for ostream
#include <vector>    // for vector, operator==

#include <hibf/detail/prefixes.hpp> // for header, merged_bin
#include <hibf/platform.hpp>

namespace hibf::layout
{

// Currently, the layout is structured by user bin.
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
            stream << prefix::header << prefix::merged_bin << '_';
            auto it = object.previous_TB_indices.begin();
            auto end = object.previous_TB_indices.end();
            // If not empty, we join with ';'
            if (it != end)
            {
                stream << *it;
                while (++it != end)
                    stream << ';' << *it;
            }
            stream << " max_bin_id:" << object.id;

            return stream;
        }
    };

    struct user_bin
    {
        std::vector<size_t> previous_TB_indices{}; // previous technical bin indices which refer to merged bin indices.
        size_t storage_TB_id{};                    // id of the technical bin that the user bin is actuallly stored in
        size_t number_of_technical_bins{};         // 1 == single bin, >1 == split_bin
        size_t idx{};                              // The index of the user bin corresponding to the order in data

        user_bin() = default;
        user_bin(user_bin const &) = default;
        user_bin(user_bin &&) = default;
        user_bin & operator=(user_bin const &) = default;
        user_bin & operator=(user_bin &&) = default;
        ~user_bin() = default;

        user_bin(size_t const idx_,
                 std::vector<size_t> const & previous_TB_indices_,
                 size_t const number_of_technical_bins_,
                 size_t const storage_TB_id_) :
            previous_TB_indices{previous_TB_indices_},
            storage_TB_id{storage_TB_id_},
            number_of_technical_bins{number_of_technical_bins_},
            idx{idx_}
        {}

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

    size_t top_level_max_bin_id{};
    std::vector<max_bin> max_bins{};
    std::vector<user_bin> user_bins{};
};

} // namespace hibf::layout
