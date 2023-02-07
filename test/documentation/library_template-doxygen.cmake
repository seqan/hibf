# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_HAVE_DOT})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (LIBRARY_TEMPLATE_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (LIBRARY_TEMPLATE_DOXYGEN_HAVE_DOT "YES")
endif ()

### Use mathjax instead of latex to render formulas.
find_package (LATEX QUIET)

if (${LATEX_FOUND})
    set (LIBRARY_TEMPLATE_DOXYGEN_USE_MATHJAX "NO")
else ()
    set (LIBRARY_TEMPLATE_DOXYGEN_USE_MATHJAX "YES")
endif ()

### Number of threads to use for dot. Doxygen's default is 0 (all threads).
set (LIBRARY_TEMPLATE_DOXYGEN_DOT_NUM_THREADS "0")

### Configure doc/developer targets.
set (LIBRARY_TEMPLATE_DOXYGEN_SOURCE_DIR "${LIBRARY_TEMPLATE_CLONE_DIR}")
set (LIBRARY_TEMPLATE_DOXYFILE_IN ${LIBRARY_TEMPLATE_DOXYGEN_INPUT_DIR}/library_template_doxygen_cfg.in)
set (LIBRARY_TEMPLATE_FOOTER_HTML_IN ${LIBRARY_TEMPLATE_DOXYGEN_INPUT_DIR}/library_template_footer.html.in)
# DoxygenLayout.xml.in is created by library_template-doxygen-layout.cmake
set (LIBRARY_TEMPLATE_LAYOUT_IN ${CMAKE_CURRENT_BINARY_DIR}/DoxygenLayout.xml.in)

option (LIBRARY_TEMPLATE_USER_DOC "Create build target and test for user documentation." ON)
option (LIBRARY_TEMPLATE_DEV_DOC "Create build target and test for developer documentation." ON)
option (LIBRARY_TEMPLATE_VERCEL_PREVIEW_DOC "Is this a preview build by vercel.com?" OFF)

if (LIBRARY_TEMPLATE_VERCEL_PREVIEW_DOC)
    set (LIBRARY_TEMPLATE_DOXYGEN_USE_MATHJAX "YES")
    set (LIBRARY_TEMPLATE_DOXYGEN_DOT_NUM_THREADS "2")
    set (LIBRARY_TEMPLATE_DOXYFILE_OPTION_POWERED_BY_VERCEL
         "HTML_EXTRA_FILES       += ${LIBRARY_TEMPLATE_DOXYGEN_SOURCE_DIR}/test/documentation/.vercel/powered-by-vercel.svg"
    )
    set (LIBRARY_TEMPLATE_FOOTER_HTML_OPTION_POWERED_BY_VERCEL
         "<li class='footer'><a href='https://vercel.com/?utm_source=seqan&utm_campaign=oss'><img class='footer' src='$relpath^powered-by-vercel.svg' height='31' alt='Powered by Vercel' style='width:unset;'/></a></li>"
    )
endif ()

### Download and extract cppreference-doxygen-web.tag.xml for std:: documentation links
set (LIBRARY_TEMPLATE_DOXYGEN_STD_TAGFILE "${PROJECT_BINARY_DIR}/cppreference-doxygen-web.tag.xml")
include (ExternalProject)
ExternalProject_Add (
    download-cppreference-doxygen-web-tag
    URL "https://github.com/PeterFeicht/cppreference-doc/releases/download/v20220201/html-book-20220201.tar.xz"
    URL_HASH SHA256=b41960e7ec9c5433b31b1b33db5854f97770ae49535c81e7647f618003102d1a
    TLS_VERIFY ON
    DOWNLOAD_DIR "${PROJECT_BINARY_DIR}"
    DOWNLOAD_NAME "html-book.tar.xz"
    DOWNLOAD_NO_EXTRACT YES
    BINARY_DIR "${PROJECT_BINARY_DIR}"
    CONFIGURE_COMMAND /bin/sh -c "xzcat html-book.tar.xz | tar -xf - cppreference-doxygen-web.tag.xml"
    BUILD_COMMAND rm "html-book.tar.xz"
    INSTALL_COMMAND "")

### TEST HELPER

# doxygen does not show any warnings (doxygen prints warnings / errors to cerr)
set (LIBRARY_TEMPLATE_TEST_DOXYGEN_FAIL_ON_WARNINGS
     "${DOXYGEN_EXECUTABLE} > doxygen.cout 2> doxygen.cerr; cat \"doxygen.cerr\"; test ! -s \"doxygen.cerr\""
     CACHE INTERNAL "The doxygen test command")

### install helper

# make sure that prefix path is /usr/local/share/doc/library_template/
if (NOT DEFINED CMAKE_SIZEOF_VOID_P)
    # we need this to suppress GNUInstallDirs AUTHOR_WARNING:
    #   CMake Warning (dev) at /usr/share/cmake-3.19/Modules/GNUInstallDirs.cmake:223 (message):
    #     Unable to determine default CMAKE_INSTALL_LIBDIR directory because no
    #     target architecture is known.  Please enable at least one language before
    #     including GNUInstallDirs.
    set (CMAKE_SIZEOF_VOID_P 8)
endif ()
include (GNUInstallDirs) # this is needed to prefix the install paths
