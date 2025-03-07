# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

if (NOT DEFINED CMAKE_CXX_STANDARD)
    set (CMAKE_CXX_STANDARD 23)
endif ()

if (NOT DEFINED CMAKE_CXX_STANDARD_REQUIRED)
    set (CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

if (NOT DEFINED CMAKE_CXX_EXTENSIONS)
    set (CMAKE_CXX_EXTENSIONS OFF)
endif ()

# LTO support.
include (CheckIPOSupported)
check_ipo_supported (
    RESULT HIBF_TEST_HAS_LTO
    OUTPUT HIBF_TEST_HAS_LTO_OUTPUT
    LANGUAGES CXX)
if (HIBF_TEST_HAS_LTO)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()
