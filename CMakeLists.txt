# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)

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

include (cmake/cxx_config.cmake)

set (HIBF_HEADER_PATH "${HIBF_SOURCE_DIR}/include")
include (cmake/configuration.cmake)
add_subdirectory (src)

option (INSTALL_HIBF "Enable installation of HIBF. (Projects embedding HIBF may want to turn this OFF.)" ON)
if (INSTALL_HIBF)
    include (GNUInstallDirs)
    include (cmake/cpack.cmake)
    include (cmake/install.cmake)
endif ()
