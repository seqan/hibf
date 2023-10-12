# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# HIBF. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.10)

# Force alignment of benchmarked loops so that numbers are reliable.
# For large loops and erratic seeming bench results the value might
# have to be adapted or the option deactivated.
option (HIBF_BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

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

# hibf::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** hibf tests
if (NOT TARGET hibf::test)
    add_library (hibf_test INTERFACE)

    # GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
            target_compile_options (hibf_test INTERFACE "-Wno-interference-size")
        endif ()
    endif ()

    target_link_libraries (hibf_test INTERFACE "seqan::hibf")
    target_include_directories (hibf_test INTERFACE "${HIBF_TEST_INCLUDE_DIR}")
    add_library (hibf::test ALIAS hibf_test)
endif ()

# hibf::test::performance specifies required flags, includes and libraries
# needed for performance test cases in hibf/test/performance
if (NOT TARGET hibf::test::performance)
    add_library (hibf_test_performance INTERFACE)
    target_link_libraries (hibf_test_performance INTERFACE "hibf::test" "benchmark_main" "benchmark")

    if (HIBF_BENCHMARK_ALIGN_LOOPS)
        target_compile_options (hibf_test_performance INTERFACE "-falign-loops=32")
    endif ()

    add_library (hibf::test::performance ALIAS hibf_test_performance)
endif ()

# hibf::test::unit specifies required flags, includes and libraries
# needed for unit test cases in hibf/test/unit
if (NOT TARGET hibf::test::unit)
    add_library (hibf_test_unit INTERFACE)
    target_link_libraries (hibf_test_unit INTERFACE "hibf::test" "gtest_main" "gtest")
    add_library (hibf::test::unit ALIAS hibf_test_unit)
endif ()

# hibf::test::header specifies required flags, includes and libraries
# needed for header test cases in hibf/test/header
if (NOT TARGET hibf::test::header)
    add_library (hibf_test_header INTERFACE)
    target_link_libraries (hibf_test_header INTERFACE "hibf::test::unit")
    target_link_libraries (hibf_test_header INTERFACE "hibf::test::performance")
    target_compile_definitions (hibf_test_header INTERFACE -DHIBF_DISABLE_DEPRECATED_WARNINGS)
    target_compile_definitions (hibf_test_header INTERFACE -DHIBF_HEADER_TEST)
    add_library (hibf::test::header ALIAS hibf_test_header)
endif ()

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in hibf.
# ----------------------------------------------------------------------------

include (hibf_test_component)
include (hibf_test_files)
include (hibf_require_benchmark)
include (hibf_require_test)
include (hibf_add_subdirectories)

get_filename_component (HIBF_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
add_subdirectory ("${HIBF_ROOT_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/hibf_lib")
target_compile_options (hibf PUBLIC "-pedantic" "-Wall" "-Wextra" "-Werror")
