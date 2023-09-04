#include <vector>    // for vector

#include <hibf/config.hpp>                        // for config
#include <hibf/layout/compute_fpr_correction.hpp> // for compute_fpr_correction
#include <hibf/layout/data_store.hpp>             // for data_store
#include <hibf/layout/execute.hpp>                // for execute
#include <hibf/layout/hierarchical_binning.hpp>   // for hierarchical_binning

namespace seqan::hibf::layout
{

size_t execute(seqan::hibf::config const & config, seqan::hibf::layout::data_store & data)
{
    data.fpr_correction =
        compute_fpr_correction({.fpr = config.maximum_false_positive_rate, // prevent clang-format
                                .hash_count = config.number_of_hash_functions,
                                .t_max = config.tmax});

    return seqan::hibf::layout::hierarchical_binning{data, config}.execute();
}

} // namespace seqan::hibf::layout
