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

find_path (HIBF_MODULE_PATH "hibf-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system/")
list (APPEND CMAKE_MODULE_PATH "${HIBF_MODULE_PATH}")

include (hibf-config-version)

if (CMAKE_VERSION VERSION_LESS 3.12)
    project (hibf
             LANGUAGES CXX
             VERSION "${HIBF_PROJECT_VERSION}")
else ()
    project (hibf
             LANGUAGES CXX
             VERSION "${HIBF_PROJECT_VERSION}"
             DESCRIPTION "HIBF -- the modern C++ library for sequence analysis" # since cmake 3.9
             HOMEPAGE_URL "https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter" # since cmake 3.12
    )
endif ()

find_package (HIBF "${HIBF_PROJECT_VERSION}" REQUIRED HINTS ${HIBF_MODULE_PATH})

option (INSTALL_HIBF "Enable installation of hibf. (Projects embedding hibf may want to turn this OFF.)" ON)

if (INSTALL_HIBF)
    include (hibf-install)
    include (hibf-package)
endif ()
