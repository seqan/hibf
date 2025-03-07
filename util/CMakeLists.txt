# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (hibf_util CXX)

# Dependency: seqan::hibf
get_filename_component (HIBF_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
include (${HIBF_ROOT_DIR}/cmake/cxx_config.cmake)
add_subdirectory ("${HIBF_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/hibf_lib")

# Dependency: Sharg
include (FetchContent)
option (SHARG_NO_TDL "" ON)
FetchContent_Declare (
    sharg
    GIT_REPOSITORY https://github.com/seqan/sharg-parser
    GIT_TAG main)
FetchContent_MakeAvailable (sharg)

add_executable (fpr_correction_check fpr_correction_check.cpp)
target_link_libraries (fpr_correction_check seqan::hibf sharg::sharg)
