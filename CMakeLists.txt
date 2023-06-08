# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# HIBF. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.4)

include (cmake/version.cmake)

project (HIBF
         LANGUAGES CXX
         VERSION "${HIBF_VERSION}"
         DESCRIPTION "HIBF -- the modern C++ library for sequence analysis" # since cmake 3.9
         HOMEPAGE_URL "https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter" # since cmake 3.12
)

if (HIBF_RELEASE_CANDIDATE_VERSION VERSION_GREATER "0")
    set (HIBF_VERSION "${HIBF_VERSION}-rc.${HIBF_RELEASE_CANDIDATE_VERSION}")
endif ()

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel." FORCE)
endif ()

include (GNUInstallDirs)
include (cmake/cpack.cmake)

set (HIBF_HEADER_PATH "${HIBF_SOURCE_DIR}/include")
include (cmake/configuration.cmake)
add_subdirectory (src)
