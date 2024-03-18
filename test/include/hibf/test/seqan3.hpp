// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Adds seqan3 to the test environment.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <hibf/platform.hpp>

// SeqAn 3 [optional]
/*!\def HIBF_HAS_SEQAN3
 * \brief Whether SeqAn3 library is available or not.
 * \ingroup test
 */
#if __has_include(<seqan3/version.h>)

#    include <seqan3/version.h>

// Require minimum version of seqan2
#    if SEQAN3_VERSION_MAJOR == 3 && SEQAN3_VERSION_MINOR >= 0
#        define HIBF_HAS_SEQAN3 1
#    endif

#endif
