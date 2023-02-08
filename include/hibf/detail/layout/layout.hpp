#pragma once

#include <string>
#include <vector>

namespace hibf
{

//!\brief Stores all information needed to construct an HIBF
struct layout
{
    std::string layout_str{};

    //!\brief Store the max bi nid for each merged bin. Merged bins are identified by there indices.
    std::vector<std::pair<std::vector<size_t>, size_t>> merged_bin_max_ids{};
};

} // namespace hibf
