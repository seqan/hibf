#pragma once

#include <iostream>
#include <set>

#include <hibf/detail/configuration.hpp>
#include <hibf/detail/layout/compute_fp_correction.hpp>
#include <hibf/detail/layout/hierarchical_binning.hpp>
#include <hibf/detail/layout/output.hpp>

namespace hibf
{

int execute(hibf::configuration & config, hibf::data_store & data)
{
   if (config.disable_estimate_union)
        config.disable_rearrangement = true;

    if (config.tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (size_t number_samples = data.kmer_counts.size();
            number_samples >= 1ULL << 32) // sqrt is bigger than uint16_t
            throw std::invalid_argument{"Too many samples. Please set a tmax (see help via `-hh`)."}; // GCOVR_EXCL_LINE
        else
            config.tmax = hibf::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(number_samples))));
    }
    else if (config.tmax % 64 != 0)
    {
        config.tmax = hibf::next_multiple_of_64(config.tmax);
        std::cerr << "[HIBF LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config.tmax << ".\n";
    }

    data.fp_correction = hibf::layout::compute_fp_correction(config.false_positive_rate, config.num_hash_functions, config.tmax);

    size_t max_hibf_id = hibf::layout::hierarchical_binning{data, config}.execute();

    return max_hibf_id;
}

} // namespace hibf
