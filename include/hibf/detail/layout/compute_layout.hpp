#pragma once

#include <hibf/config.hpp>               // for config
#include <hibf/detail/layout/layout.hpp> // for layout
#include <vector>

#include <hibf/config.hpp>
#include <hibf/detail/layout/layout.hpp>
#include <hibf/detail/sketch/hyperloglog.hpp>

namespace hibf::layout
{

layout compute_layout(config const & hibf_config,
                      std::vector<size_t> & kmer_counts,
                      std::vector<sketch::hyperloglog> & sketches);

layout compute_layout(config const & hibf_config);

} // namespace hibf::layout
