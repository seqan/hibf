#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size
#include <hibf/misc/print.hpp>               // for print, print_t

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
    seqan::hibf::print(result); // [0,0,0,1,0,0,0,0,0,0,0,0]

    // Calling `increase_bin_number_to` invalidates the agent.
    ibf.increase_bin_number_to(seqan::hibf::bin_count{60u});

    // So make sure to construct a new membership_agent.
    agent = ibf.membership_agent();
}
