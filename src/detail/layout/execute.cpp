#include <cinttypes> // for uint16_t
#include <cstddef>   // for size_t
#include <iostream>  // for char_traits, operator<<, basic_ostream, cerr
#include <stdexcept> // for invalid_argument
#include <vector>    // for vector

#include <hibf/config.hpp>                             // for config
#include <hibf/detail/data_store.hpp>                  // for data_store
#include <hibf/detail/layout/execute.hpp>              // for execute
#include <hibf/detail/layout/hierarchical_binning.hpp> // for hierarchical_binning

namespace hibf
{

size_t execute(hibf::config const & config, hibf::data_store & data)
{
    return hibf::layout::hierarchical_binning{data, config}.execute();
}

} // namespace hibf
