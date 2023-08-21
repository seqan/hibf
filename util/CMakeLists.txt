# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10...3.22)
project (hibf_util CXX)

# Dependency: seqan::hibf
get_filename_component (HIBF_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
add_subdirectory ("${HIBF_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/hibf_lib")

# Dependency: Sharg
include (FetchContent)
FetchContent_Declare (
    sharg
    URL "https://github.com/seqan/sharg-parser/releases/download/1.1.1/sharg-1.1.1-Source.tar.xz"
    URL_HASH SHA256=7330f06501718e7871e55e5fd70d0e41472cc8b34bd0e3519f8c5547510c671c)
FetchContent_Populate (sharg)
list (APPEND CMAKE_PREFIX_PATH "${sharg_SOURCE_DIR}/build_system")
find_package (sharg 1.0 REQUIRED)

add_executable (fpr_quality fpr_quality.cpp)
target_link_libraries (fpr_quality seqan::hibf sharg::sharg)
