# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.10...3.22)
project (hibf_util CXX)

# https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
    cmake_policy (SET CMP0135 NEW)
endif ()

# Dependency: seqan::hibf
get_filename_component (HIBF_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
add_subdirectory ("${HIBF_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/hibf_lib")

# Dependency: Sharg
include (FetchContent)
FetchContent_Declare (
    sharg
    URL "https://github.com/seqan/sharg-parser/releases/download/1.1.1/sharg-1.1.1-Source.tar.xz"
    URL_HASH SHA256=7330f06501718e7871e55e5fd70d0e41472cc8b34bd0e3519f8c5547510c671c
    PATCH_COMMAND curl https://github.com/seqan/sharg-parser/commit/bbbfca65333b00317b3e99a35e64f44a5c3224d6.patch |
                  patch --strip=1 --force)
FetchContent_Populate (sharg)
list (APPEND CMAKE_PREFIX_PATH "${sharg_SOURCE_DIR}/build_system")
find_package (sharg 1.0 REQUIRED)

add_executable (fpr_correction_check fpr_correction_check.cpp)
target_link_libraries (fpr_correction_check seqan::hibf sharg::sharg)
