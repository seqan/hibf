#pragma once

#include <iostream>
#include <set>

#include <seqan3/search/dream_index/detail/configuration.hpp>
#include <seqan3/search/dream_index/detail/layout/aggregate_by.hpp>
#include <seqan3/search/dream_index/detail/layout/hierarchical_binning.hpp>
#include <seqan3/search/dream_index/detail/layout/output.hpp>

namespace seqan3::hibf
{

int execute(seqan3::hibf::configuration & config, seqan3::hibf::data_store & data)
{
    if (config.rearrange_user_bins)
        config.estimate_union = true;

    if (config.tmax % 64 != 0)
    {
        config.tmax = seqan3::hibf::next_multiple_of_64(config.tmax);
        std::cerr << "[CHOPPER LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config.tmax << ".\n";
    }

    data.compute_fp_correction(config.false_positive_rate, config.num_hash_functions, config.tmax);

    seqan3::hibf::hibf_statistics global_stats{config, data.fp_correction, data.kmer_counts};
    data.stats = &global_stats.top_level_ibf;
    size_t dummy{};

    size_t max_hibf_id = seqan3::hibf::hierarchical_binning{data, config}.execute();

    if (config.output_verbose_statistics)
    {
        global_stats.print_header_to(std::cout);
        global_stats.print_summary_to(dummy, std::cout);
    }

    return max_hibf_id;
}

} // namespace seqan3::hibf
