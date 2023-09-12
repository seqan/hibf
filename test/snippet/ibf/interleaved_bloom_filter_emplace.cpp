#include <hibf/interleaved_bloom_filter.hpp>

int main()
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{12u}, seqan::hibf::bin_size{8192u}};

    // Insert the values `126`, `712` and `237` into bins `0`, `3` and `9` of the Interleaved Bloom Filter.
    ibf.emplace(126, seqan::hibf::bin_index{0u});
    ibf.emplace(712, seqan::hibf::bin_index{3u});
    ibf.emplace(237, seqan::hibf::bin_index{9u});
}
