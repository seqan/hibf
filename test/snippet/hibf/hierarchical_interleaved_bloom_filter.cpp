#include <cinttypes>  // for int64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iostream>   // for operator<<, basic_ostream, cout, char_traits
#include <vector>     // for vector

#include <hibf/config.hpp>                                // for config, insert_iterator
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/misc/print.hpp>

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

    seqan::hibf::config config{.input_fn = my_input, .number_of_user_bins = 2};

    // construct the HIBF
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

    // query the HIBF
    std::vector<size_t> query{1u, 2u, 3u};
    std::vector<size_t> query2{8u, 9u, 10u};

    auto agent = hibf.membership_agent();              // you need an agent for efficient queries
    auto & result = agent.membership_for(query, 2u);   // both user bins have hashes 1,2,3
    seqan::hibf::print(result);                        // [1,0]
    agent.sort_results();                              // Results can also be sorted
    seqan::hibf::print(result);                        // [0,1]
    auto & result2 = agent.membership_for(query2, 2u); // only user bin 0 has hashes 8,9,10
    seqan::hibf::print(result2);                       // [0]
}
