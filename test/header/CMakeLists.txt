# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (hibf_header_test CXX)

include (../hibf-test.cmake)

# We compile each header twice in separate compilation units, where the second compilation omits the header guard (check
# for cyclic includes). Each alone is sufficient to test that the header is functional, but both are needed to check for
# link errors, which can happen if the header accidentally defines a variable, e.g. a global or class static member.
# Furthermore this tests that header guards are working by including the same header twice.
#
# example invocation:
#     hibf_header_test (hibf "<path>/include/hibf" "test.hpp|/folder/|test2.hpp")
#
# \param component        The component name, will create the target `${component}_header_test`
# \param header_base_path The base path to the header files
# \param exclude_regex    A regular expression on the header file paths that excludes them from the header test. For
#                         regex syntax see https://cmake.org/cmake/help/v3.15/command/string.html#regular-expressions.
#
# \sa Modified version from Bio-Formats
# https://github.com/openmicroscopy/bioformats/blob/d3bb33eeda23e81f78fd25f658bfc14a4363805f/cpp/cmake/HeaderTest.cmake#L81-L113
#
# \sa https://en.cppreference.com/w/cpp/language/storage_duration#external_linkage
#
# ======== Detailed Explanation ========
#
# `test/header/generate_header_source.cmake` generates for each given header file two `.cpp` files.
#
# Let header be `include/hibf/search/fm_index/bi_fm_index_cursor.hpp` and component be `hibf`.
#
# * `<build_dir>/hibf_header_test/hibf/search/fm_index/bi_fm_index_cursor.hpp-header-guard.cpp`
#   which checks if the header implements a header guard.
#
#   This will be checked by including the header twice, if the header guard is missing the second include will produce
#   redeclaration errors, i.e.
#
#   ```c++
#   #include <hibf/search/fm_index/bi_fm_index_cursor.hpp>
#   #include <hibf/search/fm_index/bi_fm_index_cursor.hpp> // redeclaration errors if bi_fm_index_cursor.hpp has no header-guards
#   ```
#
# * `<build_dir>/hibf_header_test/hibf/search/fm_index/bi_fm_index_cursor.hpp-no-self-include.cpp`
#   which checks if the header will not be included itself (cyclic include). For this we copy the content of the header
#   file into the source file and if any included header includes the header itself again, we would get redeclaration
#   errors, i.e.
#   ```c++
#   // This is file bi_fm_index_cursor.hpp without #pragma once
#
#   // ...
#   #include <hibf/search/fm_index/fm_index.hpp> // redeclaration errors if fm_index.hpp includes bi_fm_index_cursor.hpp
#   ```
#
# Each `.cpp` file generates an object file which will be linked into one big binary. This ensures that we don't have
# any leaking symbols (e.g. https://en.cppreference.com/w/cpp/language/storage_duration#external_linkage).
macro (hibf_header_test component header_base_path exclude_regex)
    set (target "${component}_header_test")

    # finding all *.hpp files relative from the current directory (e.g. /test/)
    # The resulting list is normalized to `header_base_path` that means concatenating
    # "${header_base_path}/header_files[i]" will result in an absolute path to the file
    #
    # Example output:
    #   hibf/alphabet/adaptation/all.hpp
    #   hibf/alphabet/adaptation/char.hpp
    #   hibf/alphabet/adaptation/concept.hpp
    #   hibf/alphabet/adaptation/uint.hpp
    #   hibf/alphabet/all.hpp
    #   hibf/alphabet/dna5_detail.hpp <- will be filtered out
    #   ....
    hibf_test_files (header_files "${header_base_path}" "*.hpp;*.h")

    # filter out headers
    if (NOT ";${exclude_regex};" STREQUAL ";;")
        list (FILTER header_files EXCLUDE REGEX "${exclude_regex}")
    endif ()

    file (WRITE "${PROJECT_BINARY_DIR}/${target}.cpp" "")
    add_executable (${target} ${PROJECT_BINARY_DIR}/${target}.cpp)
    target_link_libraries (${target} hibf::test::header)
    add_test (NAME "header/${target}" COMMAND ${target})

    foreach (header ${header_files})
        hibf_test_component (header_test_name "${header}" TEST_NAME)
        hibf_test_component (header_target_name "${header}" TARGET_UNIQUE_NAME)

        foreach (header_sub_test "header-guard" "no-self-include")
            set (header_target_source
                 "${PROJECT_BINARY_DIR}/${target}_files/${header_test_name}.hpp-${header_sub_test}.cpp")
            set (header_target "${target}--${header_target_name}-${header_sub_test}")

            string (REPLACE "-" "__" header_test_name_safe "${target}, ${header_target}")

            # we use add_custom_command to detect changes to a header file, which will update the generated source file
            add_custom_command (OUTPUT "${header_target_source}"
                                COMMAND "${CMAKE_COMMAND}" #
                                        "-DHEADER_FILE_ABSOLUTE=${header_base_path}/${header}"
                                        "-DHEADER_FILE_INCLUDE=${header}"
                                        "-DHEADER_TARGET_SOURCE=${header_target_source}"
                                        "-DHEADER_TEST_NAME_SAFE=${header_test_name_safe}"
                                        "-DHEADER_COMPONENT=${component}" #
                                        "-DHEADER_SUB_TEST=${header_sub_test}" #
                                        "-P" "${CMAKE_CURRENT_SOURCE_DIR}/generate_header_source.cmake"
                                DEPENDS "${header_base_path}/${header}"
                                        "${CMAKE_CURRENT_SOURCE_DIR}/generate_header_source.cmake")

            add_library (${header_target} OBJECT "${header_target_source}")
            # Link hibf::test first, even though it is also linked by hibf::test:header.
            # Without this, the compile options of hibf::test are appended, e.g.,
            # `-Wno-error=... -Werror -Wall`, instead of `-Werror -Wall -Wno-error=...`
            target_link_libraries (${header_target} hibf::test hibf::test::header)
            target_sources (${target} PRIVATE $<TARGET_OBJECTS:${header_target}>)
        endforeach ()
    endforeach ()

    unset (target)
    unset (header_files)
    unset (header_test_name)
    unset (header_test_name_safe)
    unset (header_target_name)
    unset (header_target_source)
    unset (header_target)
endmacro ()

CPMGetPackage (benchmark)
if (benchmark_ADDED)
    set_target_properties (benchmark PROPERTIES CXX_STANDARD 20)
endif ()
CPMGetPackage (googletest)

# note: hibf/version.hpp is one of the only header that is not required to have a hibf/platform.hpp include
hibf_header_test (hibf "${HIBF_SOURCE_DIR}/include" "hibf/version.hpp|hibf/contrib")
hibf_header_test (hibf_test "${HIBF_SOURCE_DIR}/test/include" "hibf/test/iterator_test_template.hpp")
