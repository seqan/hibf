#include <cstddef>  // for size_t
#include <iostream> // for operator<<, basic_ostream, cout, char_traits

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, binning_bitvector, bin_index, bin_count

void print(seqan::hibf::binning_bitvector const & vector)
{
    std::cout << '[';

    if (vector.size() != 0u)
    {
        for (size_t i = 0u; i < vector.size() - 1u; ++i)
            std::cout << vector[i] << ',';
        std::cout << vector[vector.size() - 1u];
    }

    std::cout << "]\n";
}

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});

    // Query the Interleaved Bloom Filter. Note that there may be false positive results!
    // A `1` at position `i` indicates the (probable) presence of the query in bin `i`.
    // Capture the result by reference to avoid copies.
    auto agent = ibf.membership_agent();
    auto & result = agent.bulk_contains(712);
    print(result); // [0,0,0,1,0,0,0,0,0,0,0,0]

    // Calling `increase_bin_number_to` invalidates the agent.
    ibf.increase_bin_number_to(seqan::hibf::bin_count{60u});

    // So make sure to construct a new membership_agent.
    agent = ibf.membership_agent();
}
