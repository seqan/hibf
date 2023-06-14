#pragma once

#include <cstddef> // for size_t
#include <vector>  // for vector

#include <hibf/platform.hpp>

namespace hibf::layout
{

struct fpr_correction_parameters
{
    double fpr{};
    size_t hash_count{};
    size_t t_max{};
};

/*!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
 * \sa https://godbolt.org/z/zTj1v9W94
 */
std::vector<double> compute_fpr_correction(fpr_correction_parameters const & params);

} // namespace hibf::layout
