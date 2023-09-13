#include <hibf/config.hpp>                                // for insert_iterator, config
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter

void print(std::vector<int64_t> const & vector)
{
    std::cout << '[';

    if (!vector.empty())
    {
        for (size_t i = 0u; i < vector.size() - 1u; ++i)
            std::cout << vector[i] << ',';
        std::cout << vector.back();
    }

    std::cout << "]\n";
}

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
    print(result);                                     // [1,0]
    agent.sort_results();                              // Results can also be sorted
    print(result);                                     // [0,1]
    auto & result2 = agent.membership_for(query2, 2u); // only user bin 0 has hashes 8,9,10
    print(result2);                                    // [0]
}
