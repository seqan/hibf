#pragma once

#include <string_view> // for string_view

#include <hibf/platform.hpp>

namespace hibf::prefix
{
/* These prefixes are for writing the layout file

 * It is structured like this:
 *
 * [0) Possibly metadata added by chopper/raptor-layout]
 * 1) Metadata: the hibf config
 * 2) Layout header: max bin ids for the merged bins
 * 3) Layout content: Assignment of user bin idx to technical bin idx
 *
 * And marked like this:
 * [0) First character is @; Start and End of meta data should be marked accordingly to (1)]
 * 1) First character is @; Start and End are marked by @HIBF_CONFIG and @HIBF_CONFIG_END respectively
 * 2) First character is #;
 * 3) No mark, plain content.
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
 *
 * ``
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

} // namespace hibf::prefix
