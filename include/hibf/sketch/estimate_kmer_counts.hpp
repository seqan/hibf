#pragma once

#include <algorithm> // for fill_n, max
#include <cstddef>   // for size_t
#include <vector>    // for vector

#include <hibf/sketch/hyperloglog.hpp> // for hyperloglog

namespace seqan::hibf::sketch
{

inline void estimate_kmer_counts(std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                 std::vector<size_t> & kmer_counts)
{
    kmer_counts.resize(sketches.size());

    for (size_t i = 0; i < sketches.size(); ++i)
        kmer_counts[i] = sketches[i].estimate();
}

} // namespace seqan::hibf::sketch
