#pragma once

#include <string_view> // for string_view

#include <hibf/platform.hpp>

namespace seqan::hibf::prefix
{

/*!\defgroup hibf_layout_prefixes Prefixes
 * \ingroup hibf_layout
 * \brief Prefixes for writing the layout file.
 * \details
 * 1. Optional metadata added by chopper/raptor-layout
 * 2. Metadata: the hibf config
 * 3. Layout header: max bin ids for the merged bins
 * 4. Layout content: Assignment of user bin idx to technical bin idx
 *
 * And marked like this:
 * 1. First character is @; Start and End of meta data should be marked accordingly to (1)
 * 2. First character is @; Start and End are marked by \@HIBF_CONFIG and \@HIBF_CONFIG_END respectively
 * 3. First character is #;
 * 4. No mark, plain content.
 *
 * Example:
 *
 * ```
 * @CHOPPER_USER_BINS
 * @0 /path/to/file1.fa
 * @CHOPPER_USER_BINS_END
 * @CHOPPER_CONFIG
 * @0 k = 20
 * @CHOPPER_CONFIG_END
 * ```
 * \{
 */
constexpr std::string_view meta_header{"@"};

constexpr std::string_view meta_hibf_config_start{"@HIBF_CONFIG"};
static_assert(meta_hibf_config_start.starts_with(meta_header));

constexpr std::string_view meta_hibf_config_end{"@HIBF_CONFIG_END"};
static_assert(meta_hibf_config_end.starts_with(meta_header));

constexpr std::string_view layout_header{"#"};

constexpr std::string_view layout_top_level{"TOP_LEVEL_IBF"};

constexpr std::string_view layout_lower_level{"LOWER_LEVEL_IBF"};

constexpr std::string_view layout_fullest_technical_bin_idx{"fullest_technical_bin_idx:"};

constexpr std::string_view layout_first_header_line{"#TOP_LEVEL_IBF"};
static_assert(layout_first_header_line.starts_with(layout_header));
static_assert(layout_first_header_line.ends_with(layout_top_level));

constexpr std::string_view layout_column_names{"#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS"};
static_assert(layout_column_names.starts_with(layout_header));
//!\}

} // namespace seqan::hibf::prefix
