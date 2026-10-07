// SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>       // for array
#include <concepts>    // for derived_from
#include <cstddef>     // for size_t
#include <cstdint>     // for uint8_t
#include <format>      // for formatter, format_error, format_to, format
#include <functional>  // for function
#include <iosfwd>      // for ostream, istream
#include <optional>    // for optional
#include <string>      // for string
#include <string_view> // for string_view
#include <vector>      // for operator==, vector

#include <hibf/layout/prefixes.hpp> // for layout_fullest_technical_bin_idx, layout_header, layout_lower_level
#include <hibf/platform.hpp>

namespace seqan::hibf
{

struct config;

} // namespace seqan::hibf

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

    /*!\brief Returns the number of levels of the described HIBF.
     * \returns `0` if there are no user bins, `1` if there is only the top-level IBF, and so on.
     * \details
     * Only the user bins are considered.
     */
    [[nodiscard]] size_t number_of_levels() const;

    /*!\name Validation
     * \{
     */
    /*!\brief A finding of seqan::hibf::layout::layout::validate.
     * \details
     * Errors describe layouts that cannot be used to build an HIBF. Warnings describe layouts that can be used, but
     * probably do not match the intention. Notes are informational.
     */
    struct diagnostic
    {
        //!\brief How grave a finding is.
        enum class severity : uint8_t
        {
            note,
            warning,
            error
        };

        //!\brief The kind of finding. The errors are listed roughly in the order the checks are run.
        enum class code : uint8_t
        {
            // errors
            empty_layout,               //!< The layout contains no user bins.
            no_technical_bins,          //!< A user bin occupies zero technical bins.
            technical_bin_overflow,     //!< An IBF would have more than `2^64 - 64` technical bins.
            user_bin_out_of_range,      //!< A user bin index is not below `config.number_of_user_bins`.
            duplicate_user_bin,         //!< A user bin index occurs more than once.
            missing_user_bin,           //!< A user bin index is missing from the layout.
            overlapping_technical_bins, //!< Within an IBF, a technical bin is used more than once.
            top_level_max_bin_entry,    //!< `max_bins` contains an entry for the top-level IBF.
            max_bin_without_ibf,        //!< `max_bins` contains an entry for an IBF that does not exist.
            duplicate_max_bin,          //!< `max_bins` contains more than one entry for an IBF.
            invalid_max_bin,            //!< A max bin is neither a merged bin nor a user bin's first technical bin.
            max_bin_exceeds_root,       //!< A max bin spans more than `next_multiple_of_64(#TBs of the Root-IBF)` TBs.
            missing_max_bin,            //!< `max_bins` contains no entry for a lower-level IBF.
            // warnings
            technical_bin_exceeds_tmax, //!< An IBF, with its empty bins, exceeds `next_multiple_of_64(tmax)` TBs.
            single_bin_ibf,             //!< A lower-level IBF only contains a single user bin or merged bin.
            unexpected_empty_bins,      //!< An IBF does not end with the empty TBs `config.empty_bin_fraction` implies.
            // notes
            empty_technical_bins, //!< An IBF has empty technical bins between its used technical bins.
        };

        severity level{}; //!< How grave the finding is.
        code what{};      //!< The kind of finding.
        //!\brief The merged bins on the path to the affected IBF, empty for the Root-IBF. std::nullopt: no IBF.
        std::optional<std::vector<size_t>> ibf{};
        std::optional<size_t> user_bin{};      //!< The affected user bin index, if any.
        std::optional<size_t> technical_bin{}; //!< The affected technical bin of `ibf`, if any.
        std::string message{};                 //!< A human-readable description.

        //!\brief Prints `[HIBF LAYOUT <SEVERITY>] <message>`, like `std::format("{}", object)`.
        // needs a template (instead of using std::ostream directly) to be able to only include <iosfwd>
        // `diagnostic_t` is always `diagnostic`. As a template parameter, it defers the check of the format string
        // until std::formatter<diagnostic>, which can only be specialised after this class, is known.
        template <typename stream_type, std::same_as<diagnostic> diagnostic_t>
            requires std::derived_from<stream_type, std::ostream>
        friend stream_type & operator<<(stream_type & stream, diagnostic_t const & object)
        {
            stream << std::format("{}", object);
            return stream;
        }
    };

    /*!\brief Receives the findings of seqan::hibf::layout::layout::validate.
     * \details
     * The handler is called once per finding, in the order the checks run. Checking stops at the first error, so a
     * handler receives at most one error, and warnings and notes only if there is no error. A handler may throw to
     * abort the check; the layout is not modified.
     *
     * seqan::hibf::layout::layout::validate uses seqan::hibf::layout::layout::throw_on_error by default. With an empty
     * handler, it reports nothing and only returns whether the layout is valid.
     *
     * seqan::hibf::layout::layout::throw_on_error only uses seqan::hibf::layout::layout::diagnostic::level and
     * seqan::hibf::layout::layout::diagnostic::message. The other members are meant for user-defined handlers,
     * for example, to filter by seqan::hibf::layout::layout::diagnostic::code or to locate the affected IBF.
     *
     * ### Custom handlers
     *
     * Collect all findings:
     * \snippet test/snippet/hibf/layout/layout_validate.cpp collect
     *
     * Print all findings:
     * \snippet test/snippet/hibf/layout/layout_validate.cpp print
     *
     * Ignore a kind of finding:
     * \snippet test/snippet/hibf/layout/layout_validate.cpp ignore
     *
     * Treat warnings as errors:
     * \snippet test/snippet/hibf/layout/layout_validate.cpp strict
     *
     * Locate the affected IBF, user bin and technical bin:
     * \snippet test/snippet/hibf/layout/layout_validate.cpp locate
     */
    using diagnostic_handler = std::function<void(diagnostic const &)>;

    /*!\brief The default diagnostic handler of seqan::hibf::layout::layout::validate.
     * \param[in] finding The finding to handle.
     * \throws std::invalid_argument if `finding` is an error. The description is `std::format("{}", finding)`.
     * \details
     * Prints warnings to std::cerr. Prints notes to std::cerr in debug builds of the library (without `NDEBUG`).
     */
    static void throw_on_error(diagnostic const & finding);

    /*!\brief Checks whether the layout describes a consistent HIBF for the given configuration.
     * \param[in] config  The configuration the layout was or will be used with.
     * \param[in] handler Called for each finding. May be empty.
     * \returns `true` if there is no error, `false` otherwise.
     * \throws std::invalid_argument describing the first error, if `handler` is the default,
     *         seqan::hibf::layout::layout::throw_on_error.
     * \details
     * The layout must contain each user bin index in `[0, config.number_of_user_bins)` exactly once. The user bins may
     * be in any order.
     *
     * The errors are checked first, roughly in the order of seqan::hibf::layout::layout::diagnostic::code.
     * Checking stops at the first error. Warnings and notes are only reported if there is no error. They are reported
     * per IBF, so a note about one IBF may precede a warning about another IBF.
     * The build rounds the technical bins of each IBF, plus the empty bins for `config.empty_bin_fraction`, up to a
     * multiple of 64. This number is compared to `config.tmax`, unless `config.tmax` is `0`. The technical bins after
     * the last used technical bin must be exactly the empty bins `config.empty_bin_fraction` implies, i.e., none if it
     * is `0`.
     *
     * ### Example
     * \snippet test/snippet/hibf/layout/layout_validate.cpp validate
     */
    bool validate(config const & config, diagnostic_handler const & handler = throw_on_error) const;
    //!\}

    size_t top_level_max_bin_id{};
    std::vector<max_bin> max_bins{};
    std::vector<user_bin> user_bins{};

    bool operator==(layout const &) const = default;
};

} // namespace seqan::hibf::layout

//!\brief Formats a seqan::hibf::layout::layout::diagnostic as `[HIBF LAYOUT <SEVERITY>] <message>`.
template <>
struct std::formatter<seqan::hibf::layout::layout::diagnostic>
{
    constexpr formatter() = default;                              //!< Defaulted.
    constexpr formatter(formatter const &) = default;             //!< Defaulted.
    constexpr formatter & operator=(formatter const &) = default; //!< Defaulted.
    constexpr formatter(formatter &&) = default;                  //!< Defaulted.
    constexpr formatter & operator=(formatter &&) = default;      //!< Defaulted.
    constexpr ~formatter() = default;                             //!< Defaulted.

    /*!\brief Accepts only an empty format specification, e.g., `"{}"` or `"{:}"`.
     * \throws std::format_error for any other format specification, e.g., `"{:>30}"`.
     */
    template <typename parse_context_t>
    constexpr auto parse(parse_context_t & context)
    {
        auto const it = context.begin();
        if (it != context.end() && *it != '}')
            throw std::format_error{"seqan::hibf::layout::layout::diagnostic does not support format specifications."};
        return it;
    }

    /*!\brief Writes `[HIBF LAYOUT <SEVERITY>] <message>`.
     * \throws std::format_error if the severity is not valid.
     */
    template <typename format_context_t>
    auto format(seqan::hibf::layout::layout::diagnostic const & diagnostic, format_context_t & context) const
    {
        using severity = seqan::hibf::layout::layout::diagnostic::severity;
        if (diagnostic.level > severity::error)
            throw std::format_error{"Invalid seqan::hibf::layout::layout::diagnostic::severity."};

        constexpr std::array<std::string_view, 3> names{"NOTE", "WARNING", "ERROR"};
        return std::format_to(context.out(),
                              "[HIBF LAYOUT {}] {}",
                              names[static_cast<size_t>(diagnostic.level)],
                              diagnostic.message);
    }
};
