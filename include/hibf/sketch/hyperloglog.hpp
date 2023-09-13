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

#include <algorithm> // for copy
#include <array>     // for array
#include <cinttypes> // for uint64_t, uint8_t
#include <cstddef>   // for size_t
#include <iosfwd>    // for istream, ostream
#include <vector>    // for vector

#include <hibf/contrib/aligned_allocator.hpp> // for aligned_allocator

namespace seqan::hibf::sketch
{

/*!\brief HyperLogLog estimates.
 * \ingroup hibf_sketch
 * \details
 * Copied from Hideaki Ohno and adjusted/improved by Felix Droop
 * \see https://github.com/hideo55/cpp-HyperLogLog
 */
class hyperloglog
{
public:
    /*!\brief Constructor
     * \param[in] b bit width (register size will be 2 to the b power).
     *            This value must be in the range[4,30].Default value is 5.
     *
     * \throws std::invalid_argument if the argument b is out of range.
     */
    hyperloglog(uint8_t b = 5);

    /*!\brief Adds element to the estimator
     * \param[in] str string to add
     * \param[in] len length of string
     */
    void add(char const * str, uint64_t len);

    /*!\brief Estimates cardinality value.
     * \returns Estimated cardinality value.
     */
    double estimate() const;

    /*!\brief Merges the estimate from 'other' into this object.
     * \param[in] other HyperLogLog instance to be merged
     * \details
     * The number of registers in each must be the same.
     */
    void merge(hyperloglog const & other);

    /*!\brief Merges the estimate from 'other' into this object
     * \param[in] other HyperLogLog instance to be merged
     * \returns estimated cardinality of the new merged sketch.
     * \details
     * The number of registers in each must be the same.
     * This function is implemented using SIMD instructions.
     * \warning This function is undefined bevahior if this.b_ == 4
     */
    double merge_and_estimate_SIMD(hyperloglog const & other);

    /*!\brief Clears all internal registers.
     */
    void clear();

    /*!\brief Returns size of register.
     */
    uint64_t registerSize() const
    {
        return m_;
    }

    /*!\brief Exchanges the content of the instance.
     * \param[in,out] rhs Another HyperLogLog instance
     */
    void swap(hyperloglog & rhs);

    /*!\brief Dumps the current status to a stream.
     * \param[in,out] os The output stream where the data is saved to
     * \throws std::runtime_error if dumping failed.
     */
    void dump(std::ostream & os) const;

    /*!\brief Restorse the status from a stream.
     * \param[in] is The input stream where the status is saved
     * \throws std::runtime_error if restoring failed.
     */
    void restore(std::istream & is);

private:
    static constexpr std::array<float, 61> exp2_rcp = []() constexpr
    {
        std::array<float, 61> arr{};
        for (size_t i = 0; i < 61; ++i)
            arr[i] = 1.0f / static_cast<float>(1ULL << i);
        return arr;
    }();

    uint64_t mask_{};                                                                  ///< mask for the rank bits
    double alphaMM_{};                                                                 ///< alpha * m^2
    float alphaMM_float_{};                                                            ///< alpha * m^2
    uint64_t m_{};                                                                     ///< register size
    uint8_t b_{};                                                                      ///< register bit width
    std::vector<uint8_t, seqan::hibf::contrib::aligned_allocator<uint8_t, 256u>> M_{}; ///< registers
};

} // namespace seqan::hibf::sketch
