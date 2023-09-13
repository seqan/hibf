#include <algorithm> // for copy
#include <cinttypes> // for uint16_t, uint8_t
#include <concepts>  // for same_as
#include <cstddef>   // for size_t
#include <iostream>  // for operator<<, basic_ostream, cout, char_traits
#include <ranges>    // for iota_view, operator==, _Iota, iota, views
#include <vector>    // for vector

#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter, bin_index, bin_count, bin_size, count...

template <typename counter_t>
    requires (std::same_as<counter_t, uint16_t> || std::same_as<counter_t, uint8_t>)
void print(seqan::hibf::counting_vector<counter_t> const & vector)
{
    std::cout << '[';

    if (!vector.empty())
    {
        for (size_t i = 0u; i < vector.size() - 1u; ++i)
            std::cout << static_cast<uint16_t>(vector[i]) << ',';
        std::cout << static_cast<uint16_t>(vector.back());
    }

    std::cout << "]\n";
}

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{8u},
                                              seqan::hibf::bin_size{8192u},
                                              seqan::hibf::hash_function_count{2u}};

    auto sequence1 = std::views::iota(0u, 20u);
    auto sequence2 = std::views::iota(10u, 30u);
    auto sequence3 = std::views::iota(25u, 35u);

    // Insert all values of sequence1 into bin 0
    for (auto && value : sequence1)
        ibf.emplace(value, seqan::hibf::bin_index{0u});

    // Insert all values of sequence2 into bin 4
    for (auto && value : sequence2)
        ibf.emplace(value, seqan::hibf::bin_index{4u});

    // Insert all values of sequence3 into bin 7
    for (auto && value : sequence3)
        ibf.emplace(value, seqan::hibf::bin_index{7u});

    auto agent = ibf.counting_agent();

    // Count all values of sequence1 for all bins
    auto & result = agent.bulk_count(sequence1); // Bind by `&` to avoid copies!
    print(result);                               // [20,0,0,0,10,0,0,0]

    // Search for specific values
    std::vector<size_t> const values{92, 1238, 812, 81273};
    print(agent.bulk_count(values));                      // [0,0,0,0,0,0,0,0]
    print(agent.bulk_count(std::views::iota(0u, 1024u))); // [20,0,0,0,20,0,0,10]

    // The default counters are 16 bit unsigned integer.
    // An optional template parameter can be used to specify the counter type
    auto agent2 = ibf.counting_agent<uint8_t>();
    // The returned counts are now 8 bit unsigned integers.
    print(agent2.bulk_count(sequence1)); // [20,0,0,0,10,0,0,0]
}
