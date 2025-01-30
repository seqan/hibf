# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (hibf_test_snippet CXX)

include (../hibf-test.cmake)
include (../cmake/diagnostics/list_unused_snippets.cmake)

add_library (snippet_main snippet_main.cpp)
target_link_libraries (snippet_main PUBLIC hibf::test gtest)

macro (hibf_snippet test_name_prefix snippet snippet_base_path)
    hibf_test_component (snippet_target_name "${snippet}" TARGET_NAME)
    hibf_test_component (snippet_test_name "${snippet}" TEST_NAME)
    hibf_test_component (snippet_target_path "${snippet}" TARGET_PATH)
    set (target "${snippet_target_name}_snippet")

    add_executable (${target} "${snippet_base_path}/${snippet}")
    target_link_libraries (${target} PUBLIC snippet_main)
    set_target_properties (${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                "${CMAKE_CURRENT_BINARY_DIR}/${snippet_target_path}")
    collect_used_snippet (${target})

    set (snippet_test_target "${test_name_prefix}/${snippet_test_name}_snippet")
    add_test (NAME "${snippet_test_target}" COMMAND ${target})

    # disable version checker, as it interferes with comparing the snippet output
    set_tests_properties ("${snippet_test_target}" PROPERTIES ENVIRONMENT HIBF_NO_VERSION_CHECK=0)

    set (snippet_compare_test_target "${test_name_prefix}/${snippet_test_name}_snippet_cmp_output")

    add_test (NAME "${snippet_compare_test_target}"
              COMMAND ${CMAKE_COMMAND} -DTARGET_FILE=$<TARGET_FILE:${target}>
                      -DSOURCE_FILE=${snippet_base_path}/${snippet} #
                      -P "${CMAKE_CURRENT_SOURCE_DIR}/compare_snippet_output.cmake")

    unset (snippet_target_name)
    unset (snippet_test_name)
    unset (snippet_target_path)
    unset (target)
endmacro ()

macro (hibf_snippets test_name_prefix snippet_base_path)
    hibf_test_files (snippet_files "${snippet_base_path}" "*.cpp")
    foreach (snippet ${snippet_files})
        hibf_snippet ("${test_name_prefix}" "${snippet}" "${snippet_base_path}")
    endforeach ()
endmacro ()

CPMGetPackage (googletest)

hibf_snippets ("snippet" "${CMAKE_CURRENT_SOURCE_DIR}")
hibf_snippets ("doc/snippet" "${CMAKE_CURRENT_SOURCE_DIR}/../../doc")

# These are the tests that are known to have non-deterministic output.
# set_tests_properties ("snippet/.../test_name_snippet_cmp_output" PROPERTIES SKIP_RETURN_CODE 1)

list_unused_snippets ("${CMAKE_CURRENT_SOURCE_DIR}")
list_unused_snippets ("${CMAKE_CURRENT_SOURCE_DIR}/../../doc")
