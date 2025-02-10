# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)

### Find doxygen and dependency to DOT tool
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_HAVE_DOT})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (HIBF_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (HIBF_DOXYGEN_HAVE_DOT "YES")
endif ()

### Number of threads to use for dot. Doxygen's default is 0 (all threads).
set (HIBF_DOXYGEN_DOT_NUM_THREADS "0")

### Configure doc/developer targets.
set (HIBF_DOXYGEN_SOURCE_DIR "${HIBF_HEADER_PATH}/..")
set (HIBF_DOXYFILE_IN ${HIBF_DOXYGEN_INPUT_DIR}/hibf_doxygen_cfg.in)
set (HIBF_FOOTER_HTML_IN ${HIBF_DOXYGEN_INPUT_DIR}/hibf_footer.html.in)
# DoxygenLayout.xml.in is created by hibf-doxygen-layout.cmake
set (HIBF_LAYOUT_IN ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml.in)

### Download and extract cppreference-doxygen-web.tag.xml for std:: documentation links
# HIBF_DOXYGEN_STD_TAGFILE can be used to point to an existing tag file (cppreference-doxygen-web.tag.xml).
# If HIBF_DOXYGEN_STD_TAGFILE is set by the user and the file exists, it will be copied.
# If HIBF_DOXYGEN_STD_TAGFILE is not set by the user, or it is set by the user, but the file does not exist,
# the tag file will be downloaded.
set (HIBF_DEFAULT_DOXYGEN_STD_TAGFILE "${PROJECT_BINARY_DIR}/cppreference-doxygen-web.tag.xml")
set (HIBF_DOXYGEN_STD_TAGFILE
     "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}"
     CACHE STRING "Path to cppreference-doxygen-web.tag.xml")
if (NOT EXISTS "${HIBF_DOXYGEN_STD_TAGFILE}" OR HIBF_DOXYGEN_STD_TAGFILE STREQUAL "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}")
    message (STATUS "Tag file will be fetched.")
    # Reset path in case it was set from the outside, but does not exist.
    set (HIBF_DOXYGEN_STD_TAGFILE "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}")
    include (ExternalProject)
    # When updating, check whether warnings in HIBF_TEST_DOXYGEN_FAIL_ON_WARNINGS are gone when removing sed filter.
    ExternalProject_Add (
        download-cppreference-doxygen-web-tag
        URL "https://github.com/PeterFeicht/cppreference-doc/releases/download/v20230810/html-book-20230810.tar.xz"
        URL_HASH SHA256=31c08e4d99e86c7f63f324d3ff5304eff2030131c4a0ac0d1e3c19c62c8ed684
        TLS_VERIFY ON
        DOWNLOAD_DIR "${PROJECT_BINARY_DIR}"
        DOWNLOAD_NAME "html-book.tar.xz"
        DOWNLOAD_NO_EXTRACT YES
        BINARY_DIR "${PROJECT_BINARY_DIR}"
        BUILD_BYPRODUCTS "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}"
        CONFIGURE_COMMAND /bin/sh -c "xzcat html-book.tar.xz | tar -xf - cppreference-doxygen-web.tag.xml"
        BUILD_COMMAND rm "html-book.tar.xz"
        INSTALL_COMMAND "")
else ()
    message (STATUS "Copying existing tag file: ${HIBF_DOXYGEN_STD_TAGFILE}")
    # Copy tag file such that it is present in the built documentation. This is useful if the documentation is
    # subsequently deployed to a server.
    add_custom_target (download-cppreference-doxygen-web-tag)
    add_custom_command (TARGET download-cppreference-doxygen-web-tag
                        POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy "${HIBF_DOXYGEN_STD_TAGFILE}"
                                "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}"
                        BYPRODUCTS "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}")
    set (HIBF_DOXYGEN_STD_TAGFILE "${HIBF_DEFAULT_DOXYGEN_STD_TAGFILE}")
endif ()

### TEST HELPER

# doxygen does not show any warnings (doxygen prints warnings / errors to cerr)
set (HIBF_TEST_DOXYGEN_FAIL_ON_WARNINGS
     "${DOXYGEN_EXECUTABLE} > doxygen.cout 2> doxygen.cerr; cat \"doxygen.cerr\"; test ! -s \"doxygen.cerr\""
     CACHE INTERNAL "The doxygen test command")

### install helper

# make sure that prefix path is /usr/local/share/doc/hibf/
if (NOT DEFINED CMAKE_SIZEOF_VOID_P)
    # we need this to suppress GNUInstallDirs AUTHOR_WARNING:
    #   CMake Warning (dev) at /usr/share/cmake-3.19/Modules/GNUInstallDirs.cmake:223 (message):
    #     Unable to determine default CMAKE_INSTALL_LIBDIR directory because no
    #     target architecture is known.  Please enable at least one language before
    #     including GNUInstallDirs.
    set (CMAKE_SIZEOF_VOID_P 8)
endif ()
include (GNUInstallDirs) # this is needed to prefix the install paths
