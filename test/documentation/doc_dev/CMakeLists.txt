# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

message (STATUS "Configuring devel doc.")

set (LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (LIBRARY_TEMPLATE_DOXYGEN_EXCLUDE_SYMBOLS "")
set (LIBRARY_TEMPLATE_DOXYGEN_PREDEFINED_NDEBUG "")
set (LIBRARY_TEMPLATE_DOXYGEN_ENABLED_SECTIONS "DEV")
set (LIBRARY_TEMPLATE_SHOW_DEV_DOCS "YES")

configure_file (${LIBRARY_TEMPLATE_DOXYFILE_IN} ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file (${LIBRARY_TEMPLATE_FOOTER_HTML_IN} ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/library_template_footer.html)
configure_file (${LIBRARY_TEMPLATE_LAYOUT_IN} ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/DoxygenLayout.xml)

add_custom_target (doc_dev ALL
                   COMMAND ${DOXYGEN_EXECUTABLE}
                   WORKING_DIRECTORY ${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}
                   BYPRODUCTS html/
                   DEPENDS download-cppreference-doxygen-web-tag
                   COMMENT "Generating developer API documentation with Doxygen"
                   VERBATIM)

# Install doc_dev documentation in ./install_doc_dev folder
# cmake --install . --prefix install_doc_dev --component doc_dev
install (DIRECTORY "${LIBRARY_TEMPLATE_DOXYGEN_OUTPUT_DIR}/html"
         DESTINATION "${CMAKE_INSTALL_DOCDIR}"
         COMPONENT doc-dev
         EXCLUDE_FROM_ALL
         PATTERN "*.md5" EXCLUDE
         PATTERN "*.map" EXCLUDE
         PATTERN "formula.repository" EXCLUDE)

### Enable testing

enable_testing ()

message (STATUS "Add dev doc test.")
add_test (NAME doc_dev_no_warnings_test
          COMMAND bash -c "${LIBRARY_TEMPLATE_TEST_DOXYGEN_FAIL_ON_WARNINGS}"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
