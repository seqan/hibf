// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

/*!\file
 * \brief Provides version macros and global variables.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

//!\brief The major version as MACRO.
#define HIBF_VERSION_MAJOR 1
//!\brief The minor version as MACRO.
#define HIBF_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define HIBF_VERSION_PATCH 0
//!\brief The release candidate number. 0 means stable release, >= 1 means release candidate.
#define HIBF_RELEASE_CANDIDATE 1

//!\brief The full version as MACRO (number).
#define HIBF_VERSION                                                                                       \
    (HIBF_VERSION_MAJOR * 10000 + HIBF_VERSION_MINOR * 100 + HIBF_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define HIBF_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define HIBF_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH)                                              \
    HIBF_VERSION_CSTRING_HELPER_STR(MAJOR)                                                                 \
    "." HIBF_VERSION_CSTRING_HELPER_STR(MINOR) "." HIBF_VERSION_CSTRING_HELPER_STR(PATCH)

#if (HIBF_RELEASE_CANDIDATE > 0)
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define HIBF_RELEASE_CANDIDATE_HELPER(RC) "-rc." HIBF_VERSION_CSTRING_HELPER_STR(RC)
#else
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define HIBF_RELEASE_CANDIDATE_HELPER(RC) ""
#endif

//!\brief The full version as null terminated string.
#define HIBF_VERSION_CSTRING                                                                               \
    HIBF_VERSION_CSTRING_HELPER_FUNC(HIBF_VERSION_MAJOR,                                       \
                                                 HIBF_VERSION_MINOR,                                       \
                                                 HIBF_VERSION_PATCH)                                       \
    HIBF_RELEASE_CANDIDATE_HELPER(HIBF_RELEASE_CANDIDATE)

namespace hibf
{

//!\brief The major version.
constexpr uint8_t hibf_version_major = HIBF_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t hibf_version_minor = HIBF_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t hibf_version_patch = HIBF_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t hibf_version = HIBF_VERSION;

//!\brief The full version as null terminated string.
constexpr char const * hibf_version_cstring = HIBF_VERSION_CSTRING;

} // namespace hibf

#undef HIBF_VERSION_CSTRING_HELPER_STR
#undef HIBF_VERSION_CSTRING_HELPER_FUNC
#undef HIBF_RELEASE_CANDIDATE_HELPER
