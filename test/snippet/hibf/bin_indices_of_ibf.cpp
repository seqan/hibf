#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

// For this example we have two input fasta files with the following three sequences (= user bins):
// example1.fasta:
// ```fasta
// >chr1
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
// ```
// example2.fasta:
// ```fasta
// >chr2
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTCATTAA
// >chr3
// AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
// ```

// 1-level-IBF: |chr1|chr2,chr3|
// 2-level-IBF:      |chr2|chr3|

template <std::ranges::range range_t>
void print_range(range_t && range)
{
    std::cerr << '[';
    auto range_begin = std::ranges::begin(range);
    auto range_end = std::ranges::end(range);
    auto it = range_begin;
    for (; it != std::ranges::prev(range_end, 1, range_begin); ++it)
        std::cerr << *it << ',';
    if (range_begin != range_end)
        std::cerr << *it;
    std::cerr << ']';
}

int main()
{
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{};

    hibf.user_bins.set_ibf_count(2);

    hibf.user_bins.bin_indices_of_ibf(0) = {1, 2, 3};
    hibf.user_bins.bin_indices_of_ibf(1) = {2, 3};

    std::cerr << "User bin indices of 1-level-IBF: ";
    print_range(hibf.user_bins.bin_indices_of_ibf(0));
    std::cerr << '\n';

    std::cerr << "User bin indices of 2-level-IBF: ";
    print_range(hibf.user_bins.bin_indices_of_ibf(1));
    std::cerr << '\n';
}

// Prints:
// User bin indices of 1-level-IBF: [1,2,3]
// User bin indices of 2-level-IBF: [2,3]
