# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.12...3.30)

function (generate_include_dependencies_impl)
    cmake_parse_arguments (
        "" #
        "" #
        "TARGET;TARGET_INTERNAL_DEPENDENCY_MAKE_FILE;HIBF_INCLUDE_DIR;TARGET_DEPENDENCIES_FILE" #
        "TARGET_CYCLIC_DEPENDING_INCLUDES" #
        ${ARGN})

    if (NOT EXISTS "${_TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}")
        return ()
    endif ()

    # File content may look like this:
    #
    # alphabet/nucleotide/CMakeFiles/dna4_test.dir/depend.make
    #
    # ```
    # # CMAKE generated file: DO NOT EDIT!
    # # Generated by "Unix Makefiles" Generator, CMake Version 3.20
    #
    # utility/views/CMakeFiles/zip_test.dir/zip_test.cpp.o: \
    #  /hibf/include/hibf/core/platform.hpp \
    #  /hibf/include/hibf/std/algorithm \
    # ```

    # read in file and filter out linebreaks
    file (STRINGS "${_TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}" header_files)

    # store original file for diagnostics
    set (header_files_raw "${header_files}")

    # filter out "\;" as they would escape semicolons which are the separators for cmake list elements
    string (REPLACE "\;" ";" header_files "${header_files}")

    # only use lines that contain a hibf include
    list (FILTER header_files INCLUDE REGEX "${_HIBF_INCLUDE_DIR}/hibf")

    # filter out object files, e.g., discard "utility/views/CMakeFiles/zip_test.dir/zip_test.cpp.o: " in the line
    # `utility/views/CMakeFiles/zip_test.dir/zip_test.cpp.o: /hibf/include/hibf/core/platform.hpp`
    list (TRANSFORM header_files REPLACE "^.+: " "")

    if (NOT header_files)
        # The pre-processing step that generates the dependency file did not produce it.
        # We will warn about this.
        #
        # There might be the following reasons:
        #
        # 1. The generation of the dependency file changed (this happened once with cmake 3.20)
        # 2. The header does not contain any hibf include.
        message (AUTHOR_WARNING "no hibf header files contained: ${_TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}\n"
                                "file: ${header_files_raw}")
    endif ()
    unset (header_files_raw)

    set (script "")
    foreach (header_file ${header_files})
        # remove leading and trailing whitespaces
        string (STRIP "${header_file}" header_file)

        file (RELATIVE_PATH relative_header_file "${_HIBF_INCLUDE_DIR}/hibf" "${header_file}")
        get_include_target (include_depended_target HEADER "${relative_header_file}")

        string (APPEND script "add_include_target (\"${include_depended_target}\")\n")
        # exclude known cyclic header
        if (NOT "${include_depended_target}" IN_LIST _TARGET_CYCLIC_DEPENDING_INCLUDES)
            string (APPEND script "add_dependencies (${_TARGET} \"${include_depended_target}\")\n")
        else ()
            string (APPEND script "# add_dependencies (${_TARGET} \"${include_depended_target}\") # known cycle\n")
        endif ()
        string (APPEND script "\n")
    endforeach ()

    file (WRITE "${_TARGET_DEPENDENCIES_FILE}" "${script}")
endfunction ()

if (CMAKE_SCRIPT_MODE_FILE)
    list (APPEND CMAKE_MODULE_PATH "${HIBF_TEST_CMAKE_MODULE_DIR}")
    include (include_dependencies/add_include_dependencies)
    include (hibf_test_component)

    message (STATUS "Generate include dependencies of target ${TARGET}")

    generate_include_dependencies_impl (
        # e.g. dna4_test
        TARGET "${TARGET}"
        # e.g. alphabet/nucleotide/CMakeFiles/dna4_test.dir/depend.make
        TARGET_INTERNAL_DEPENDENCY_MAKE_FILE
            "${TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}"
            # e.g. /hibf-repo/include
            HIBF_INCLUDE_DIR "${HIBF_INCLUDE_DIR}"
        # e.g. alphabet/nucleotide/dna4_test_dependencies.cmake (will be generated)
        TARGET_DEPENDENCIES_FILE "${TARGET_DEPENDENCIES_FILE}"
        TARGET_CYCLIC_DEPENDING_INCLUDES "${TARGET_CYCLIC_DEPENDING_INCLUDES}")
endif ()
