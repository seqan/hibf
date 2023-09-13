#include <cinttypes> // for uint16_t
#include <cstddef>   // for size_t
#include <iostream>  // for operator<<, basic_ostream, cout, char_traits
#include <vector>    // for vector

#include <hibf/interleaved_bloom_filter.hpp> // for counting_vector, interleaved_bloom_filter, bin_index, bin_count

void print(seqan::hibf::counting_vector<uint16_t> const & vector)
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
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(126, seqan::hibf::bin_index{3u});
    ibf.emplace(126, seqan::hibf::bin_index{9u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});

    // The counting_vector must be at least as big as the number of bins.
    seqan::hibf::counting_vector<uint16_t> counts(12, 0);

    auto agent = ibf.membership_agent();

    counts += agent.bulk_contains(712); // `counts` contains the number of occurrences of 712 in each bin.
    print(counts);                      // prints [0,0,0,1,0,0,0,0,0,0,0,0]

    counts += agent.bulk_contains(237); // `counts` contains the number of occurrences of 712 and 237 in each bin.
    print(counts);                      // prints [0,0,0,1,0,0,0,0,0,1,0,0]

    counts += agent.bulk_contains(126); // `counts` contains the number of occurrences of 712, 237 and 126 in each bin.
    print(counts);                      // prints [1,0,0,2,0,0,0,0,0,2,0,0]

    counts += counts; // multiple counts can also be added together
    print(counts);    // prints [2,0,0,4,0,0,0,0,0,4,0,0]
}
