#pragma once

#include <iosfwd>
#include <vector>

#include <hibf/detail/configuration.hpp>
#include <hibf/detail/layout/layout.hpp>

namespace hibf::layout
{

void write_config_to(configuration const &, std::ostream &);

void write_layout_header_to(layout const &, size_t const, std::ostream &);

void write_user_bin_line_to(layout::user_bin const &, std::vector<std::string> const &, std::ostream &);

void write_layout_content_to(layout const &, std::vector<std::string> const &, std::ostream &);

} // namespace hibf::layout
