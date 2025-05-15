# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (hibf_test_unit CXX)

include (../hibf-test.cmake)
include (GoogleTest OPTIONAL)

include (diagnostics/list_missing_unit_tests)
include (diagnostics/list_unused_unit_tests)
include (include_dependencies/add_include_dependencies)

option (HIBF_USE_INCLUDE_DEPENDENCIES
        "Build tests in an hierarchical order (by an include graph, i.e. tests with less dependencies are build first)"
        OFF)

macro (hibf_test unit_test_cpp)
    cmake_parse_arguments (HIBF_TEST "" "SDE" "CYCLIC_DEPENDING_INCLUDES" ${ARGN})

    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    hibf_test_component (target "${unit_test}" TARGET_NAME)
    hibf_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} hibf::test::unit)
    add_include_dependencies (${target} "${HIBF_TEST_CYCLIC_DEPENDING_INCLUDES}")
    collect_used_test (${target})

    if (HIBF_TEST_SDE)
        set (TARGET_WITH_SDE "${HIBF_TEST_SDE}" "-cnl" "-check_app_mt" "--" "./${target}")
        add_test (NAME "${test_name}_sde" COMMAND ${TARGET_WITH_SDE})
    else ()
        add_test (NAME "${test_name}" COMMAND ${target})
    endif ()

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

CPMGetPackage (googletest)

add_subdirectories ()

list_missing_unit_tests ()
list_unused_unit_tests ()
