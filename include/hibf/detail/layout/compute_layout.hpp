#pragma once

#include <hibf/config.hpp>               // for config
#include <hibf/detail/layout/layout.hpp> // for layout

namespace hibf::layout
{

layout compute_layout(config const & hibf_config);

} // namespace hibf::layout
