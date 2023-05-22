// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
