# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# HIBF. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.10)

# require HIBF package
find_package (HIBF REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../build_system)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

option (HIBF_TEST_BUILD_OFFLINE "Skip the update step of external projects." OFF)

# Force alignment of benchmarked loops so that numbers are reliable.
# For large loops and erratic seeming bench results the value might
# have to be adapted or the option deactivated.
option (HIBF_BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

# ----------------------------------------------------------------------------
# Custom Build types
# ----------------------------------------------------------------------------

# -DCMAKE_BUILD_TYPE=FEDORA; our library did not compile for fedora quite a few times, because of that we created this
# custom build type to emulate their flag set-up.
# We omitted:
#   -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1
#   -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1
#   set -mtune=native
#   and -fcf-protection=check
# See https://src.fedoraproject.org/rpms/redhat-rpm-config/blob/rawhide/f/buildflags.md for an overview
set (CMAKE_CXX_FLAGS_FEDORA
     "-O2 -flto -ffat-lto-objects -fexceptions -g -grecord-gcc-switches -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fstack-protector-strong -m64 -mtune=native -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection=check"
)

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
    target_compile_options (hibf_test INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror")

    # GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
            target_compile_options (hibf_test INTERFACE "-Wno-interference-size")
        endif ()
    endif ()

    target_link_libraries (hibf_test INTERFACE "hibf::hibf" "pthread")
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

# hibf::test::coverage specifies required flags, includes and libraries
# needed for coverage test cases in hibf/test/coverage
if (NOT TARGET hibf::test::coverage)
    add_library (hibf_test_coverage INTERFACE)
    target_compile_options (hibf_test_coverage INTERFACE "--coverage" "-fprofile-arcs" "-ftest-coverage")
    # -fprofile-abs-path requires at least gcc8, it forces gcov to report absolute instead of relative paths.
    # gcovr has trouble detecting the headers otherwise.
    # ccache is not aware of this option, so it needs to be skipped with `--ccache-skip`.
    find_program (CCACHE_PROGRAM ccache)
    if (CCACHE_PROGRAM)
        target_compile_options (hibf_test_coverage INTERFACE "--ccache-skip" "-fprofile-abs-path")
    else ()
        target_compile_options (hibf_test_coverage INTERFACE "-fprofile-abs-path")
    endif ()
    target_link_libraries (hibf_test_coverage INTERFACE "hibf::test::unit" "gcov")
    add_library (hibf::test::coverage ALIAS hibf_test_coverage)
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
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (HIBF_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND HIBF_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND HIBF_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND HIBF_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND HIBF_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in hibf.
# ----------------------------------------------------------------------------

include (hibf_test_component)
include (hibf_test_files)
include (hibf_require_ccache)
include (hibf_require_benchmark)
include (hibf_require_test)
include (hibf_add_subdirectories)
