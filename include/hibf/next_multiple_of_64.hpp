#pragma once

#include <cstddef> // for size_t

#include <hibf/platform.hpp>

namespace seqan::hibf
{

/*!\brief Returns the smallest integer that is greater or equal to `value` and a multiple of 64. Returns 0 for value 0.
 * \param[in] value The Input value that is smaller or equal to the return value.
 */
[[nodiscard]] constexpr size_t next_multiple_of_64(size_t const value) noexcept
{
    return ((value + 63) >> 6) << 6;
}

} // namespace seqan::hibf
