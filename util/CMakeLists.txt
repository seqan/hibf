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
CPMAddPackage (URI "gh:seqan/sharg-parser#main"
               SYSTEM TRUE
               OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
                       "SHARG_NO_TDL ON")

add_executable (fpr_correction_check fpr_correction_check.cpp)
target_link_libraries (fpr_correction_check seqan::hibf sharg::sharg)

add_executable (hash_collisions hash_collisions.cpp)
target_link_libraries (hash_collisions seqan::hibf sharg::sharg)

set (CMAKE_REQUIRED_QUIET 1)
check_cxx_source_runs ("#include <generator>\nint main(){}" HIBF_HAS_COROUTINE_GENERATOR)
if (HIBF_HAS_COROUTINE_GENERATOR)
    # Dependency: fmt
    CPMAddPackage (URI "gh:fmtlib/fmt#11.2.0"
                   SYSTEM TRUE
                   OPTIONS "FMT_INSTALL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING")

    add_executable (inspect inspect/inspect.cpp)
    target_link_libraries (inspect seqan::hibf sharg::sharg fmt::fmt)
else ()
    # Available starting with GCC 14.
    # Clang's libc++ does not support std::generator yet (Clang-20).
    message (STATUS " Not building inspect because std::generator is not available.")
endif ()
