# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file provides functionality common to the different test modules used by
# HIBF. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.20...3.31)

# Force alignment of benchmarked loops so that numbers are reliable.
# For large loops and erratic seeming bench results the value might
# have to be adapted or the option deactivated.
option (HIBF_BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

include (CheckIPOSupported)
check_ipo_supported (
    RESULT HIBF_TEST_HAS_LTO
    OUTPUT HIBF_TEST_HAS_LTO_OUTPUT
    LANGUAGES CXX)
if (HIBF_TEST_HAS_LTO)
    set (CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

get_filename_component (HIBF_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

option (HIBF_POST_INSTALL_TEST "Tests should use installed library." OFF)
if (HIBF_POST_INSTALL_TEST)
    find_package (hibf CONFIG REQUIRED)
else ()
    add_subdirectory ("${HIBF_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/hibf_lib")
    target_compile_options (hibf PUBLIC "-pedantic" "-Wall" "-Wextra" "-Werror")

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        # Warn about failed return value optimization.
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
            target_compile_options (hibf PUBLIC "-Wnrvo")
        endif ()
    endif ()
endif ()

set (CPM_INDENT "  CMake Package Manager CPM: ")
include (${HIBF_ROOT_DIR}/cmake/CPM.cmake)
CPMUsePackageLock (${HIBF_ROOT_DIR}/cmake/package-lock.cmake)

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

find_path (HIBF_TEST_INCLUDE_DIR
           NAMES hibf/test/tmp_directory.hpp
           HINTS "${CMAKE_CURRENT_LIST_DIR}/include/")
find_path (HIBF_TEST_CMAKE_MODULE_DIR
           NAMES hibf_test_component.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${HIBF_TEST_CMAKE_MODULE_DIR}")

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in hibf.
# ----------------------------------------------------------------------------

enable_testing ()

# hibf::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** hibf tests
if (NOT TARGET hibf::test)
    add_library (hibf_test INTERFACE)

    target_link_libraries (hibf_test INTERFACE "seqan::hibf")
    target_include_directories (hibf_test INTERFACE "${HIBF_TEST_INCLUDE_DIR}")
    add_library (hibf::test ALIAS hibf_test)
endif ()

# hibf::test::performance specifies required flags, includes and libraries
# needed for performance test cases in hibf/test/performance
if (NOT TARGET hibf::test::performance)
    add_library (hibf_test_performance INTERFACE)
    target_link_libraries (hibf_test_performance INTERFACE "hibf::test" "benchmark::benchmark_main")

    if (HIBF_BENCHMARK_ALIGN_LOOPS)
        target_compile_options (hibf_test_performance INTERFACE "-falign-loops=32")
    endif ()

    add_library (hibf::test::performance ALIAS hibf_test_performance)
endif ()

# hibf::test::unit specifies required flags, includes and libraries
# needed for unit test cases in hibf/test/unit
if (NOT TARGET hibf::test::unit)
    add_library (hibf_test_unit INTERFACE)
    target_link_libraries (hibf_test_unit INTERFACE "hibf::test" "GTest::gtest_main")
    add_library (hibf::test::unit ALIAS hibf_test_unit)
endif ()

# hibf::test::header specifies required flags, includes and libraries
# needed for header test cases in hibf/test/header
if (NOT TARGET hibf::test::header)
    add_library (hibf_test_header INTERFACE)
    target_link_libraries (hibf_test_header INTERFACE "hibf::test::unit")
    target_link_libraries (hibf_test_header INTERFACE "hibf::test::performance")
    target_compile_options (hibf_test_header INTERFACE "-Wno-unused-function" "-Wno-unused-const-variable")
    target_compile_definitions (hibf_test_header INTERFACE -DHIBF_DISABLE_DEPRECATED_WARNINGS)
    target_compile_definitions (hibf_test_header INTERFACE -DHIBF_HEADER_TEST)
    add_library (hibf::test::header ALIAS hibf_test_header)
endif ()

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in hibf.
# ----------------------------------------------------------------------------

include (hibf_test_component)
include (hibf_test_files)
include (hibf_add_subdirectories)
