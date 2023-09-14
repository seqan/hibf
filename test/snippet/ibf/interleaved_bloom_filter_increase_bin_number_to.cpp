#include <cstddef>  // for size_t
#include <iostream> // for operator<<, basic_ostream, cout, char_traits

#include <hibf/interleaved_bloom_filter.hpp> // for binning_bitvector, interleaved_bloom_filter, bin_index, bin_count
#include <hibf/misc/print.hpp>

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});

    ibf.increase_bin_number_to(seqan::hibf::bin_count{18u});
    // Be sure to get the agent after `increase_bin_number_to` as it invalidates all agents!
    auto agent = ibf.membership_agent();

    // The content of the bins which were already present before the resize does not change
    seqan::hibf::print(agent.bulk_contains(126)); // [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    seqan::hibf::print(agent.bulk_contains(712)); // [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    seqan::hibf::print(agent.bulk_contains(237)); // [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
}
