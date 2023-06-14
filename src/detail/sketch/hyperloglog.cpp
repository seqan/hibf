#include <array>     // for array
#include <cassert>   // for assert
#include <cinttypes> // for uint8_t, uint32_t
#include <cmath>     // for log
#include <cstddef>   // for size_t
#include <vector>    // for vector

#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator
#include <hibf/detail/sketch/hyperloglog.hpp> // for hyperloglog

#include <x86/avx.h>  // for simde_mm256_add_ps, simde_mm256_set_ps, simde__m256i, simde_mm...
#include <x86/avx2.h> // for simde_mm256_max_epu8

namespace hibf::sketch
{

double hyperloglog::estimate() const
{
    // compute indicator formula
    double sum = 0.0;
    for (uint8_t c : M_)
    {
        sum += exp2_rcp[c];
    }
    double estimate = alphaMM_ / sum;

    // use linear counting of zeros for small values
    if (estimate <= 2.5 * m_)
    {
        uint32_t zeros = 0;

        for (size_t i = 0; i < m_; ++i)
        {
            if (!M_[i])
                ++zeros;
        }

        if (zeros != 0u)
        {
            estimate = m_ * std::log(static_cast<double>(m_) / static_cast<double>(zeros));
        }
    }
    return estimate;
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
        uint32_t zeros = 0u;

        for (size_t i = 0; i < m_; ++i)
        {
            if (!M_[i])
                ++zeros;
        }

        if (zeros != 0u)
        {
            estimate = m_ * std::log(static_cast<double>(m_) / static_cast<double>(zeros));
        }
    }

    return estimate;
}

} // namespace hibf::sketch
