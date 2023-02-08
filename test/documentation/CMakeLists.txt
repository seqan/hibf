# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.7)

### Find hibf
include (../../build_system/hibf-config-version.cmake)
set (HIBF_VERSION "${PACKAGE_VERSION}")

project (hibf
         LANGUAGES NONE
         VERSION "${HIBF_PROJECT_VERSION}")

if (NOT EXISTS "${HIBF_INCLUDE_DIR}/hibf/version.hpp")
    message (FATAL_ERROR "Could not find HIBF. Not building documentation.")
endif ()
message (STATUS "Found: ${HIBF_INCLUDE_DIR} (${HIBF_VERSION})")

set (HIBF_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (hibf-doxygen.cmake)

include (hibf-doxygen-layout.cmake)

enable_testing ()

if (HIBF_USER_DOC)
    add_subdirectory (doc_usr)
endif ()

if (HIBF_DEV_DOC)
    add_subdirectory (doc_dev)
endif ()

include (hibf-doxygen-package.cmake)
