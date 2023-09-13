# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)

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
add_test (NAME docs_no_warnings_test
          COMMAND bash -c "${HIBF_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

include (hibf-doxygen-package.cmake)
