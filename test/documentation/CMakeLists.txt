# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.7)

include (../../cmake/version.cmake)

project (HIBF
         LANGUAGES NONE
         VERSION "${HIBF_VERSION}")

if (HIBF_RELEASE_CANDIDATE_VERSION VERSION_GREATER "0")
    set (HIBF_VERSION "${HIBF_VERSION}-rc.${HIBF_RELEASE_CANDIDATE_VERSION}")
endif ()

find_path (HIBF_HEADER_PATH
           NAMES hibf/version.hpp
           HINTS "${HIBF_SOURCE_DIR}/../../include")

if (HIBF_HEADER_PATH)
    message (STATUS "Found HIBF: ${HIBF_HEADER_PATH} (found version \"${HIBF_VERSION}\")")
else ()
    message (FATAL_ERROR "Could not find HIBF. Not building documentation.")
endif ()

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
