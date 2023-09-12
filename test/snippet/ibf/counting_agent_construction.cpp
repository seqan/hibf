#include <hibf/interleaved_bloom_filter.hpp>

int main()
{
    // Construct an Interleaved Bloom Filter to be used with the counting_agent.
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{43u},
                                              seqan::hibf::bin_size{8192u},
                                              seqan::hibf::hash_function_count{3}};

    // The counting_agent can now be constructed by calling `counting_agent` on the Interleaved Bloom Filter.
    auto agent = ibf.counting_agent();

    // Calling `increase_bin_number_to` invalidates the agent.
    ibf.increase_bin_number_to(seqan::hibf::bin_count{60u});

    // So make sure to construct a new counting_agent.
    agent = ibf.counting_agent();
}
