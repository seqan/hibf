#include <algorithm> // for fill, max
#include <array>     // for array
#include <bit>       // for countl_zero
#include <cassert>   // for assert
#include <cinttypes> // for uint8_t, uint64_t, uint32_t
#include <cmath>     // for log
#include <cstddef>   // for size_t
#include <iostream>  // for istream, ostream
#include <stdexcept> // for runtime_error, invalid_argument
#include <string>    // for operator+, to_string
#include <utility>   // for swap
#include <vector>    // for vector

#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/contrib/xxhash/xxhash.h>       // for XXH_INLINE_XXH3_64bits, XXH3_64bits
#include <hibf/sketch/hyperloglog.hpp>        // for hyperloglog

#include <x86/avx.h>  // for simde_mm256_add_ps, simde_mm256_set_ps, simde__m256i, simde_mm...
#include <x86/avx2.h> // for simde_mm256_max_epu8

namespace seqan::hibf::sketch
{

hyperloglog::hyperloglog(uint8_t b) : m_(1 << b), b_(b), M_(m_, 0)
{
    if (b < 4 || 32 < b)
        throw std::invalid_argument("bit width must be in the range [4,32] and it is " + std::to_string(b));

    M_.shrink_to_fit();
    double alpha;

    switch (m_)
    {
    case 16:
        alpha = 0.673;
        break;
    case 32:
        alpha = 0.697;
        break;
    case 64:
        alpha = 0.709;
        break;
    default:
        alpha = 0.7213 / (1.0 + 1.079 / m_);
        break;
    }

    alphaMM_ = alpha * m_ * m_;
    alphaMM_float_ = static_cast<float>(alphaMM_);
    // 64 bits where the last b are ones and the rest zeroes
    mask_ = (1 << b) - 1;
}

void hyperloglog::add(char const * str, uint64_t len)
{
    uint64_t const hash = XXH3_64bits(str, len);
    // the first b_ bits are used to distribute the leading zero counts along M_
    uint64_t const index = hash >> (64 - b_);
    // the bitwise-or with mask_ assures that we get at most 64 - b_ as value.
    // Otherwise the count for hash = 0 would be 64
    uint8_t const rank = std::countl_zero((hash << b_) | mask_) + 1;
    M_[index] = std::max(rank, M_[index]);
}

double hyperloglog::estimate() const
{
    // compute indicator formula
    double sum = 0.0;
    for (uint8_t c : M_)
        sum += exp2_rcp[c];
    double estimate = alphaMM_ / sum;

    // use linear counting of zeros for small values
    if (estimate <= 2.5 * m_)
    {
        uint32_t zeros{};

        for (size_t i = 0; i < m_; ++i)
            zeros += (M_[i] == 0u);

        if (zeros != 0u)
            estimate = m_ * std::log(static_cast<double>(m_) / static_cast<double>(zeros));
    }
    return estimate;
}

void hyperloglog::merge(hyperloglog const & other)
{
    assert(m_ == other.m_);

    for (size_t i = 0; i < m_; ++i)
    {
        if (M_[i] < other.M_[i])
        {
            M_[i] = other.M_[i];
        }
    }
}

double hyperloglog::merge_and_estimate_SIMD(hyperloglog const & other)
{
    assert(m_ == other.m_);
    assert(b_ >= 5);

    // this is safe when b_ is at least 5. Then, M_'s size in bits is
    // 2^x * 2^5 * 8 = 2^x * 256 >= 256, where x is an integer >= 1
    // also, M_ is 256 bit aligned in memory
    simde__m256i * it = reinterpret_cast<simde__m256i *>(&*(M_.begin()));
    simde__m256i const * other_it = reinterpret_cast<simde__m256i const *>(&*(other.M_.begin()));
    simde__m256i * end = reinterpret_cast<simde__m256i *>(&*(M_.end()));

    simde__m256 packed_sum = simde_mm256_set1_ps(0.0f);

    for (; it != end; ++it, ++other_it)
    {
        // this merges the registers by computing the byte-wise maximum
        *it = simde_mm256_max_epu8(*it, *other_it);

        // get pointer to iterate over the single merged registers
        uint8_t * reg_it = reinterpret_cast<uint8_t *>(it);

        // get floats with two to the power of minus the value in the merged registers and sum up
        packed_sum = simde_mm256_add_ps(packed_sum,
                                        simde_mm256_set_ps(exp2_rcp[*reg_it],
                                                           exp2_rcp[*(reg_it + 1)],
                                                           exp2_rcp[*(reg_it + 2)],
                                                           exp2_rcp[*(reg_it + 3)],
                                                           exp2_rcp[*(reg_it + 4)],
                                                           exp2_rcp[*(reg_it + 5)],
                                                           exp2_rcp[*(reg_it + 6)],
                                                           exp2_rcp[*(reg_it + 7)]));

        // repeat 3 times...
        packed_sum = simde_mm256_add_ps(packed_sum,
                                        simde_mm256_set_ps(exp2_rcp[*(reg_it + 8)],
                                                           exp2_rcp[*(reg_it + 9)],
                                                           exp2_rcp[*(reg_it + 10)],
                                                           exp2_rcp[*(reg_it + 11)],
                                                           exp2_rcp[*(reg_it + 12)],
                                                           exp2_rcp[*(reg_it + 13)],
                                                           exp2_rcp[*(reg_it + 14)],
                                                           exp2_rcp[*(reg_it + 15)]));

        packed_sum = simde_mm256_add_ps(packed_sum,
                                        simde_mm256_set_ps(exp2_rcp[*(reg_it + 16)],
                                                           exp2_rcp[*(reg_it + 17)],
                                                           exp2_rcp[*(reg_it + 18)],
                                                           exp2_rcp[*(reg_it + 19)],
                                                           exp2_rcp[*(reg_it + 20)],
                                                           exp2_rcp[*(reg_it + 21)],
                                                           exp2_rcp[*(reg_it + 22)],
                                                           exp2_rcp[*(reg_it + 23)]));

        packed_sum = simde_mm256_add_ps(packed_sum,
                                        simde_mm256_set_ps(exp2_rcp[*(reg_it + 24)],
                                                           exp2_rcp[*(reg_it + 25)],
                                                           exp2_rcp[*(reg_it + 26)],
                                                           exp2_rcp[*(reg_it + 27)],
                                                           exp2_rcp[*(reg_it + 28)],
                                                           exp2_rcp[*(reg_it + 29)],
                                                           exp2_rcp[*(reg_it + 30)],
                                                           exp2_rcp[*(reg_it + 31)]));
    }

    // sum up the 4 values in the packed SSE variable
    float sum = 0.0;
    float * sum_it = reinterpret_cast<float *>(&packed_sum);
    sum += *sum_it;
    sum += *(sum_it + 1);
    sum += *(sum_it + 2);
    sum += *(sum_it + 3);
    sum += *(sum_it + 4);
    sum += *(sum_it + 5);
    sum += *(sum_it + 6);
    sum += *(sum_it + 7);

    // compute first estimate
    double estimate = alphaMM_float_ / sum;

    // use linear counting of zeros for small values
    if (estimate <= 2.5 * m_)
    {
        uint32_t zeros{};

        for (size_t i = 0; i < m_; ++i)
            zeros += (M_[i] == 0u);

        if (zeros != 0u)
            estimate = m_ * std::log(static_cast<double>(m_) / static_cast<double>(zeros));
    }

    return estimate;
}

void hyperloglog::clear()
{
    std::fill(M_.begin(), M_.end(), 0);
}

void hyperloglog::swap(hyperloglog & rhs)
{
    std::swap(mask_, rhs.mask_);
    std::swap(alphaMM_, rhs.alphaMM_);
    std::swap(alphaMM_float_, rhs.alphaMM_float_);
    std::swap(m_, rhs.m_);
    std::swap(b_, rhs.b_);
    M_.swap(rhs.M_);
}

void hyperloglog::dump(std::ostream & os) const
{
    os.write((char *)&b_, sizeof(b_));
    os.write((char *)&M_[0], sizeof(M_[0]) * M_.size());
    os.flush();
    if (os.fail())
    {
        throw std::runtime_error("Failed to dump a HyperLogLog sketch to a file.");
    }
}

void hyperloglog::restore(std::istream & is)
{
    try
    {
        uint8_t b{};
        is.read((char *)&b, sizeof(b));
        hyperloglog tempHLL(b);
        is.read((char *)&(tempHLL.M_[0]), sizeof(M_[0]) * tempHLL.m_);
        if (is.fail())
        {
            throw std::runtime_error("Failed to restore a HyperLogLog sketch from a file.");
        }
        swap(tempHLL);
    }
    catch (std::invalid_argument const & err)
    {
        // turn the invalid argument error to a runtime error, because it is dependent on the file contents here
        throw std::runtime_error("Failed to restore a HyperLogLog sketch from a file.");
    }
}

} // namespace seqan::hibf::sketch
