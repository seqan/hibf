#pragma once

/**
 * @file hyperloglog.hpp
 * @brief HyperLogLog cardinality estimator
 * @date Created 2013/3/20, Adjusted 2021/01
 * @author Hideaki Ohno
 *
 * Copied from Hideaki Ohno (https://github.com/hideo55/cpp-HyperLogLog) and adjusted/improved by Felix Droop
 * Modified a lot for a bugfix, improvements and functional changes (64 bit hashes)
 */

#include <array>
#include <cassert>
#include <iostream>
#include <vector>

#include <hibf/contrib/aligned_allocator.hpp>
#include <hibf/contrib/xxhash/xxhash.h>
#include <hibf/platform.hpp>

namespace hibf::sketch
{

/** @class hyperloglog
 *  @brief Implement of 'HyperLogLog' estimate cardinality algorithm
 *
 * Copied from Hideaki Ohno (https://github.com/hideo55/cpp-HyperLogLog) and adjusted/improved by Felix Droop
 */
class hyperloglog
{
public:
    /**
     * Constructor
     *
     * @param[in] b bit width (register size will be 2 to the b power).
     *            This value must be in the range[4,30].Default value is 5.
     *
     * @exception std::invalid_argument the argument b is out of range.
     */
    hyperloglog(uint8_t b = 5) : m_(1 << b), b_(b), M_(m_, 0)
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

    /**
     * Adds element to the estimator
     *
     * @param[in] str string to add
     * @param[in] len length of string
     */
    void add(char const * str, uint64_t len)
    {
        uint64_t hash = XXH3_64bits(str, len);
        // the first b_ bits are used to distribute the leading zero counts along M_
        uint64_t index = hash >> (64 - b_);
        // WARNING: __builtin_clzl() only works with g++ and clang
        // the bitwise-or with mask_ assures that we get at most 64 - b_ as value.
        // Otherwise the count for hash = 0 would be 64
        uint8_t rank = __builtin_clzl((hash << b_) | mask_) + 1;
        M_[index] = std::max(rank, M_[index]);
    }

    /**
     * Estimates cardinality value.
     *
     * @return Estimated cardinality value.
     */
    double estimate() const;

    /**
     * Merges the estimate from 'other' into this object
     * The number of registers in each must be the same.
     *
     * @param[in] other HyperLogLog instance to be merged
     */
    void merge(hyperloglog const & other)
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

    /**
     * Merges the estimate from 'other' into this object
     * The number of registers in each must be the same.
     * This function is implemented using SIMD instructions.
     * WARNING: This function is undefined bevahior if this.b_ == 4
     *
     * @param[in] other HyperLogLog instance to be merged
     *
     * @return estimated cardinality of the new merged sketch
     */
    double merge_and_estimate_SIMD(hyperloglog const & other);

    /**
     * Clears all internal registers.
     */
    void clear()
    {
        std::fill(M_.begin(), M_.end(), 0);
    }

    /**
     * Returns size of register.
     *
     * @return Register size
     */
    uint64_t registerSize() const
    {
        return m_;
    }

    /**
     * Exchanges the content of the instance
     *
     * @param[in,out] rhs Another HyperLogLog instance
     */
    void swap(hyperloglog & rhs)
    {
        std::swap(mask_, rhs.mask_);
        std::swap(alphaMM_, rhs.alphaMM_);
        std::swap(alphaMM_float_, rhs.alphaMM_float_);
        std::swap(m_, rhs.m_);
        std::swap(b_, rhs.b_);
        M_.swap(rhs.M_);
    }

    /**
     * Dump the current status to a stream
     *
     * @param[out] os The output stream where the data is saved
     *
     * @exception std::runtime_error When failed to dump.
     */
    void dump(std::ostream & os) const
    {
        os.write((char *)&b_, sizeof(b_));
        os.write((char *)&M_[0], sizeof(M_[0]) * M_.size());
        os.flush();
        if (os.fail())
        {
            throw std::runtime_error("Failed to dump a HyperLogLog sketch to a file.");
        }
    }

    /**
     * Restore the status from a stream
     *
     * @param[in] is The input stream where the status is saved
     *
     * @exception std::runtime_error When failed to restore.
     */
    void restore(std::istream & is)
    {
        try
        {
            uint8_t b = 0;
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

private:
    static constexpr std::array<float, 61> exp2_rcp = []() constexpr
    {
        std::array<float, 61> arr{};
        for (size_t i = 0; i < 61; ++i)
            arr[i] = 1.0f / static_cast<float>(1ULL << i);
        return arr;
    }();

    uint64_t mask_{};                                                           ///< mask for the rank bits
    double alphaMM_{};                                                          ///< alpha * m^2
    float alphaMM_float_{};                                                     ///< alpha * m^2
    uint64_t m_{};                                                              ///< register size
    uint8_t b_{};                                                               ///< register bit width
    std::vector<uint8_t, hibf::contrib::aligned_allocator<uint8_t, 256u>> M_{}; ///< registers
};

} // namespace hibf::sketch
