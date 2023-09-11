# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.5)

include (cmake/version.cmake)

project (HIBF
         LANGUAGES CXX
         VERSION "${HIBF_VERSION}"
         DESCRIPTION "HIBF -- the modern C++ library for sequence analysis"
         HOMEPAGE_URL "https://github.com/seqan/hibf")

if (HIBF_RELEASE_CANDIDATE_VERSION VERSION_GREATER "0")
    set (HIBF_VERSION "${HIBF_VERSION}-rc.${HIBF_RELEASE_CANDIDATE_VERSION}")
endif ()

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

set (HIBF_HEADER_PATH "${HIBF_SOURCE_DIR}/include")
include (cmake/configuration.cmake)
add_subdirectory (src)

option (INSTALL_HIBF "Enable installation of HIBF. (Projects embedding HIBF may want to turn this OFF.)" ON)
if (INSTALL_HIBF)
    include (GNUInstallDirs)
    include (cmake/cpack.cmake)
    include (cmake/install.cmake)
endif ()
