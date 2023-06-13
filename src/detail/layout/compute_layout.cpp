#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iosfwd>     // for stringstream
#include <iterator>   // for inserter
#include <string>     // for operator+, to_string, basic_string, string
#include <vector>     // for vector

#include <hibf/contrib/robin_hood.hpp>   // for unordered_flat_set
#include <hibf/detail/configuration.hpp> // for configuration
#include <hibf/detail/data_store.hpp>    // for data_store
#include <hibf/detail/layout/compute_layout.hpp>
#include <hibf/detail/layout/execute.hpp>              // for execute
#include <hibf/detail/layout/layout.hpp>               // for layout
#include <hibf/detail/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/detail/sketch/hyperloglog.hpp>          // for hyperloglog

namespace hibf::layout
{

layout compute_layout(config const & hibf_config)
{
    layout resulting_layout{};

    hibf::configuration chopper_config{.sketch_bits = hibf_config.sketch_bits,
                                       .disable_sketch_output = true,
                                       .tmax = hibf_config.tmax,
                                       .num_hash_functions = hibf_config.number_of_hash_functions,
                                       .false_positive_rate = hibf_config.maximum_false_positive_rate,
                                       .alpha = hibf_config.alpha,
                                       .max_rearrangement_ratio = hibf_config.max_rearrangement_ratio,
                                       .threads = hibf_config.threads,
                                       .disable_estimate_union = hibf_config.disable_estimate_union,
                                       .disable_rearrangement = hibf_config.disable_rearrangement};

    // The output streams facilitate writing the layout file in hierarchical structure.
    // hibf::execute currently writes the filled buffers to the output file.
    std::stringstream output_buffer;
    std::stringstream header_buffer;

    std::vector<std::string> filenames{};
    std::vector<size_t> kmer_counts{};
    std::vector<sketch::hyperloglog> sketches{};

    // dummy init filenames
    filenames.resize(hibf_config.number_of_user_bins);
    for (size_t i = 0; i < hibf_config.number_of_user_bins; ++i)
        filenames[i] = "UB_" + std::to_string(i);

    // compute sketches
    sketches.resize(hibf_config.number_of_user_bins);
    kmer_counts.resize(hibf_config.number_of_user_bins);

    robin_hood::unordered_flat_set<uint64_t> kmers;
#pragma omp parallel for schedule(static) num_threads(hibf_config.threads) private(kmers)
    for (size_t i = 0; i < hibf_config.number_of_user_bins; ++i)
    {
        hibf::sketch::hyperloglog sketch(hibf_config.sketch_bits);

        kmers.clear();
        hibf_config.input_fn(i, std::inserter(kmers, kmers.begin()));

        for (auto k_hash : kmers)
            sketch.add(reinterpret_cast<char *>(&k_hash), sizeof(k_hash));

        // #pragma omp critical
        sketches[i] = sketch;
    }

    sketch::estimate_kmer_counts(sketches, kmer_counts);

    data_store store{.false_positive_rate = chopper_config.false_positive_rate,
                     .hibf_layout = &resulting_layout,
                     .kmer_counts = kmer_counts,
                     .sketches = sketches};

    size_t const max_hibf_id = hibf::execute(chopper_config, store);
    store.hibf_layout->top_level_max_bin_id = max_hibf_id;

    return *store.hibf_layout; // return layout as string for now, containing the file
}

} // namespace hibf::layout
