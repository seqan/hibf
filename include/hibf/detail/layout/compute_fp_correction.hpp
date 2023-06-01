#pragma once

#include <vector>

#include <hibf/platform.hpp>

namespace hibf::layout
{

/*!\brief Precompute f_h factors that adjust the split bin size to prevent FPR inflation due to multiple testing.
    * \sa https://godbolt.org/z/zTj1v9W94
    */
std::vector<double> compute_fp_correction(double const, size_t const, size_t const);

} // namespace hibf::layout
