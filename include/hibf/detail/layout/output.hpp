#pragma once

#include <cstddef> // for size_t
#include <iosfwd>  // for ostream
#include <string>  // for string
#include <vector>  // for vector

#include <hibf/config.hpp>               // for config
#include <hibf/detail/layout/layout.hpp> // for layout

namespace hibf::layout
{

void write_config_to(config const &, std::ostream &);

void write_layout_header_to(layout const &, size_t const, std::ostream &);

void write_user_bin_line_to(layout::user_bin const &, std::vector<std::string> const &, std::ostream &);

void write_layout_content_to(layout const &, std::vector<std::string> const &, std::ostream &);

} // namespace hibf::layout
