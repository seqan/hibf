# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library-template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

message (STATUS "Configuring user doc.")

set (LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (LIBRARY_TEMPLATE_DOXYGEN_EXCLUDE_SYMBOLS "detail library_template::simd") #/""
set (LIBRARY_TEMPLATE_DOXYGEN_PREDEFINED_NDEBUG "-NDEBUG") #/""
set (LIBRARY_TEMPLATE_DOXYGEN_ENABLED_SECTIONS "") #/"DEV"
set (LIBRARY_TEMPLATE_SHOW_DEV_DOCS "NO") #/"YES"

configure_file (${LIBRARY_TEMPLATE_DOXYFILE_IN} ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${LIBRARY_TEMPLATE_FOOTER_HTML_IN} ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/library_template_footer.html)
configure_file (${LIBRARY_TEMPLATE_LAYOUT_IN} ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/DoxygenLayout.xml)

add_custom_target (doc_usr ALL
                   COMMAND ${DOXYGEN_EXECUTABLE}
                   WORKING_DIRECTORY ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}
                   BYPRODUCTS html/
                   DEPENDS download-cppreference-doxygen-web-tag
                   COMMENT "Generating user API documentation with Doxygen"
                   VERBATIM)

# Install doc_usr documentation in ./export folder
# make DESTDIR=export install
install (DIRECTORY "${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/html"
         DESTINATION "${CMAKE_INSTALL_DOCDIR}"
         COMPONENT doc
         PATTERN "*.md5" EXCLUDE
         PATTERN "*.map" EXCLUDE
         PATTERN "formula.repository" EXCLUDE)

### Enable testing

enable_testing ()

message (STATUS "Add user doc test.")
add_test (NAME doc_usr_no_warnings_test COMMAND bash -c "${LIBRARY_TEMPLATE_TEST_DOXYGEN_FAIL_ON_WARNINGS}")
