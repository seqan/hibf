#include <cstddef>    // for size_t
#include <functional> // for function
#include <vector>     // for vector

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

int main()
{
    // 2 user bins:
    std::vector<std::vector<size_t>> hashes{{1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u}, {1u, 2u, 3u, 4u, 5u}};

    // input just passes hashes:
    auto my_input = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        for (auto const hash : hashes[user_bin_id])
            it = hash;
    };

    seqan::hibf::config config{.input_fn = my_input,     // required
                               .number_of_user_bins = 2, // required
                               .number_of_hash_functions = 2,
                               .maximum_false_positive_rate = 0.05, // recommended to adapt
                               .threads = 1,                        // recommended to adapt
                               .sketch_bits = 12,
                               .tmax = 0, // triggers default copmutation
                               .alpha = 1.2,
                               .max_rearrangement_ratio = 0.5,
                               .disable_estimate_union = false,
                               .disable_rearrangement = false};

    // construct the HIBF
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};
}
