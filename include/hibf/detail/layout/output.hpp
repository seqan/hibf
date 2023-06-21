#pragma once

#include <cstddef> // for size_t
#include <iosfwd>  // for ostream
#include <string>  // for string
#include <vector>  // for vector

#include <hibf/detail/configuration.hpp> // for configuration
#include <hibf/detail/layout/layout.hpp> // for layout

namespace hibf::layout
{

void write_config_to(configuration const &, std::ostream &);

void write_layout_header_to(layout const &, size_t const, std::ostream &);

void write_user_bin_line_to(layout::user_bin const &, std::vector<std::string> const &, std::ostream &);

void write_layout_content_to(layout const &, std::vector<std::string> const &, std::ostream &);

void write_layout_file(layout const & hibf_layout,
                       std::vector<std::string> const & filenames,
                       configuration const & config);

} // namespace hibf::layout
