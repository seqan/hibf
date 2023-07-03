#include <cinttypes> // for uint16_t
#include <cmath>     // for ceil, sqrt
#include <cstddef>   // for size_t
#include <iostream>  // for char_traits, operator<<, basic_ostream, cerr
#include <stdexcept> // for invalid_argument
#include <vector>    // for vector

#include <hibf/config.hpp>                               // for config
#include <hibf/detail/data_store.hpp>                    // for data_store
#include <hibf/detail/layout/compute_fpr_correction.hpp> // for compute_fpr_correction
#include <hibf/detail/layout/execute.hpp>                // for execute
#include <hibf/detail/layout/hierarchical_binning.hpp>   // for hierarchical_binning
#include <hibf/next_multiple_of_64.hpp>                  // for next_multiple_of_64

namespace hibf
{

size_t execute(hibf::config const & config, hibf::data_store & data)
{
    data.fpr_correction =
        layout::compute_fpr_correction({.fpr = config.maximum_false_positive_rate, // prevent clang-format
                                        .hash_count = config.number_of_hash_functions,
                                        .t_max = config.tmax});

    return hibf::layout::hierarchical_binning{data, config}.execute();
}

} // namespace hibf
