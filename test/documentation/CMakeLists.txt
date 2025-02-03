# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)

include (../../cmake/version.cmake)

project (HIBF
         LANGUAGES NONE
         VERSION "${HIBF_VERSION}")

if (HIBF_RELEASE_CANDIDATE_VERSION VERSION_GREATER "0")
    set (HIBF_VERSION "${HIBF_VERSION}-rc.${HIBF_RELEASE_CANDIDATE_VERSION}")
endif ()

find_path (HIBF_HEADER_PATH
           NAMES hibf/version.hpp
           HINTS "${HIBF_SOURCE_DIR}/../../include")

if (HIBF_HEADER_PATH)
    message (STATUS "Found HIBF: ${HIBF_HEADER_PATH} (found version \"${HIBF_VERSION}\")")
else ()
    message (FATAL_ERROR "Could not find HIBF. Not building documentation.")
endif ()

set (CPM_INDENT "CMake Package Manager CPM: ")
include (${HIBF_SOURCE_DIR}/../../cmake/CPM.cmake)
CPMUsePackageLock (${HIBF_SOURCE_DIR}/../../cmake/package-lock.cmake)

# Sets doxygen_awesome_SOURCE_DIR
CPMGetPackage (doxygen_awesome)

set (HIBF_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (hibf-doxygen.cmake)
include (hibf-doxygen-layout.cmake)

# Add documentation
message (STATUS "Configuring docs.")

set (HIBF_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (HIBF_DOXYGEN_EXCLUDE_SYMBOLS "")
set (HIBF_DOXYGEN_PREDEFINED_NDEBUG "")

configure_file (${HIBF_DOXYFILE_IN} ${HIBF_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${HIBF_FOOTER_HTML_IN} ${HIBF_DOXYGEN_OUTPUT_DIR}/hibf_footer.html)
configure_file (${HIBF_LAYOUT_IN} ${HIBF_DOXYGEN_OUTPUT_DIR}/DoxygenLayout.xml)

add_custom_target (docs ALL
                   COMMAND ${DOXYGEN_EXECUTABLE}
                   WORKING_DIRECTORY ${HIBF_DOXYGEN_OUTPUT_DIR}
                   BYPRODUCTS html/
                   DEPENDS download-cppreference-doxygen-web-tag
                   COMMENT "Generating developer API documentation with Doxygen"
                   VERBATIM)

# Install docs documentation in ./install_docs folder
# cmake --install . --prefix install_docs --component docs
install (DIRECTORY "${HIBF_DOXYGEN_OUTPUT_DIR}/html"
         DESTINATION "${CMAKE_INSTALL_DOCDIR}"
         COMPONENT docs
         EXCLUDE_FROM_ALL
         PATTERN "*.md5" EXCLUDE
         PATTERN "*.map" EXCLUDE
         PATTERN "formula.repository" EXCLUDE)

### Tests
enable_testing ()

message (STATUS "Adding docs test.")
add_test (NAME cppreference-doxygen-web-tag COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target
                                                    download-cppreference-doxygen-web-tag)
add_test (NAME docs_no_warnings_test
          COMMAND bash -c "${HIBF_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties (docs_no_warnings_test PROPERTIES DEPENDS cppreference-doxygen-web-tag)

include (hibf-doxygen-package.cmake)
