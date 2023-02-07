// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/library-template/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

/*!\file
 * \brief Provides version macros and global variables.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

//!\brief The major version as MACRO.
#define LIBRARY_TEMPLATE_VERSION_MAJOR 1
//!\brief The minor version as MACRO.
#define LIBRARY_TEMPLATE_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define LIBRARY_TEMPLATE_VERSION_PATCH 0
//!\brief The release candidate number. 0 means stable release, >= 1 means release candidate.
#define LIBRARY_TEMPLATE_RELEASE_CANDIDATE 1

//!\brief The full version as MACRO (number).
#define LIBRARY_TEMPLATE_VERSION                                                                                       \
    (LIBRARY_TEMPLATE_VERSION_MAJOR * 10000 + LIBRARY_TEMPLATE_VERSION_MINOR * 100 + LIBRARY_TEMPLATE_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH)                                              \
    LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_STR(MAJOR)                                                                 \
    "." LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_STR(MINOR) "." LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_STR(PATCH)

#if (LIBRARY_TEMPLATE_RELEASE_CANDIDATE > 0)
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define LIBRARY_TEMPLATE_RELEASE_CANDIDATE_HELPER(RC) "-rc." LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_STR(RC)
#else
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define LIBRARY_TEMPLATE_RELEASE_CANDIDATE_HELPER(RC) ""
#endif

//!\brief The full version as null terminated string.
#define LIBRARY_TEMPLATE_VERSION_CSTRING                                                                               \
    LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_FUNC(LIBRARY_TEMPLATE_VERSION_MAJOR,                                       \
                                                 LIBRARY_TEMPLATE_VERSION_MINOR,                                       \
                                                 LIBRARY_TEMPLATE_VERSION_PATCH)                                       \
    LIBRARY_TEMPLATE_RELEASE_CANDIDATE_HELPER(LIBRARY_TEMPLATE_RELEASE_CANDIDATE)

namespace library_template
{

//!\brief The major version.
constexpr uint8_t library_template_version_major = LIBRARY_TEMPLATE_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t library_template_version_minor = LIBRARY_TEMPLATE_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t library_template_version_patch = LIBRARY_TEMPLATE_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t library_template_version = LIBRARY_TEMPLATE_VERSION;

//!\brief The full version as null terminated string.
constexpr char const * library_template_version_cstring = LIBRARY_TEMPLATE_VERSION_CSTRING;

} // namespace library_template

#undef LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_STR
#undef LIBRARY_TEMPLATE_VERSION_CSTRING_HELPER_FUNC
#undef LIBRARY_TEMPLATE_RELEASE_CANDIDATE_HELPER
