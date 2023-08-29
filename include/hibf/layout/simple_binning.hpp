#pragma once

#include <cassert>   // for assert
#include <cstddef>   // for size_t
#include <stdexcept> // for runtime_error
#include <utility>   // for addressof
#include <vector>    // for vector

#include <hibf/layout/data_store.hpp>   // for data_store
#include <hibf/next_multiple_of_64.hpp> // for next_multiple_of_64

namespace seqan::hibf::layout
{

/*!\brief Distributes x Technical Bins across y User Bins while minimizing the maximal Technical Bin size
 *
 * # Terminology
 *
 * ## Technical Bin
 * \copydetails simple_binning::num_technical_bins
 *
 * ## User Bin
 * \copydetails simple_binning::num_user_bins
 *
 * ## Notation
 *  | Name    | Description                                                                 |
 *  |---------|---------------------------------------------------------------------------- |
 *  | **x**   | Number of Technical Bins (TB)                                               |
 *  | **y**   | Number of User Bins (UB)                                                    |
 *  | **b_i** | The bin size (kmer content) of Technical Bin \f$i\f$                            |
 *  | **c_j** | The kmer content of User Bin \f$j\f$                                            |
 *  | **M**   | A DP matrix that tracks the maximum technical bin size \f$\max_{i} (b_i)\f$.|
 *
 *  \attention The number of technical bins **x** must be greater that the number of user bins **y** for this algorithm.
     *            If you want to use less technical bins than user bins, see the hierarchical_binning algorithm.
 *
 * # Algorithm
 *
 * Since the size of the IBF depends on the maximal Technical Bin size, we want to minimize \f$ max_{i} (b_i)\f$.
 *
 * Let \f$r = x - y\f$ be the surplus of TBs
 *
 * ## Initialization
 * <img src="execute_init.png" align="left"  width="10%"/>
 *
 * <br><br>
 * \f$\qquad \forall_{i \in [0,r]} \quad M_{i,0} = \frac{c_0}{i + 1}\f$
 * <br><br><br><br><br><br><br><br><br><br>
 *
 * ## Recursion
 * <img src="execute_recursion.png" align="left"  width="10%"/>
 *
 * <br><br>
 * \f$\forall_{i,j} \quad M_{i,j} = \min_{i' \in [i - r - 1, i - 1]} \max(M_{i',j-1}, \frac{c_j}{i - i'})\f$
 * <br><br><br><br><br><br><br><br><br><br>
 *
 * ## Backtracking
 *
 * Assume we filled a trace matrix **T** during the computation of **M**.
 *
 * We now want to recover the number of bins **n_j** for each User Bin **j**.
 *
 * <img src="execute_backtracking.png" align="left"  width="10%"/>
 *
 * Backtracking pseudo code:
 * ```
 * // Start at the bottom-right cell.
 * j = y - 1;
 * i = x - 1;
 * n = array(y); // array of length y
 *
 * while (j > 0)
 * {
 *     next_i = T[i][j];
 *     n[j] = i - next_i;
 *     i = next_i;
 *     j = j - 1;
 * }
 * ```
 */
class simple_binning
{
private:
    //!\brief Stores all data that is needed to compute the layout, e.g. the counts, sketches and the layout::layout.
    data_store * data{nullptr};

    /*!\brief The number of User bins.
     *
     * The user may impose a structure on his sequence data in the form of *logical groups* (e.g. species).
     * When querying the IBF, the user is interested in an answer that differentiates between these groups.
     */
    size_t num_user_bins{};

    /*!\brief The number of Technical bins.
     *
     * A *Technical Bin* represents an actual bin in the binning directory.
     * In the IBF, it stores its kmers in a **single Bloom Filter** (which is interleaved with all the other BFs).
     */
    size_t num_technical_bins{};

public:
    simple_binning() = default;                                  //!< Defaulted.
    simple_binning(simple_binning const &) = delete;             //!< Deleted. Would modify same data.
    simple_binning & operator=(simple_binning const &) = delete; //!< Deleted. Would modify same data.
    simple_binning(simple_binning &&) = default;                 //!< Defaulted.
    simple_binning & operator=(simple_binning &&) = default;     //!< Defaulted.
    ~simple_binning() = default;                                 //!< Defaulted.

    /*!\brief The constructor from user bin names, their kmer counts and a configuration.
     * \param[in] data_ Stores all data that is needed to compute the layout.
     * \param[in] num_bins (optional) The number of technical bins.
     *
     * If the `num_bins` parameter is omitted or set to 0, then number of technical bins used in this algorithm
     * is automatically set to the next multiple of 64 given the number of user bins (e.g. \#UB = 88 -> \#TB = 124).
     *
     * \attention The number of technical bins must be greater or equal to the number of user bins!
     *            If you want to use less technical bins than user bins, see the hierarchical_binning algorithm.
     */
    simple_binning(data_store & data_, size_t const num_bins = 0) :
        data{std::addressof(data_)},
        num_user_bins{data->positions.size()},
        num_technical_bins{num_bins ? num_bins : next_multiple_of_64(num_user_bins)}
    {
        assert(data != nullptr);

        if (num_user_bins > num_technical_bins)
        {
            throw std::runtime_error{"You cannot have less technical bins than user bins for this simple binning "
                                     "algorithm. Please see the hierarchical_binning algorithm or increase the number "
                                     "of technical bins."};
        }
    }

    size_t get_num_technical_bins() const
    {
        return num_technical_bins;
    }

    //!\brief Executes the simple binning algorithm and layouts user bins into technical bins.
    size_t execute();
};

} // namespace seqan::hibf::layout
