# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.7)

### Find library_template
include (../../build_system/library_template-config-version.cmake)
set (LIBRARY_TEMPLATE_VERSION "${PACKAGE_VERSION}")

project (library_template
         LANGUAGES NONE
         VERSION "${LIBRARY_TEMPLATE_PROJECT_VERSION}")

if (NOT EXISTS "${LIBRARY_TEMPLATE_INCLUDE_DIR}/library_template/version.hpp")
    message (FATAL_ERROR "Could not find Library_Template. Not building documentation.")
endif ()
message (STATUS "Found: ${LIBRARY_TEMPLATE_INCLUDE_DIR} (${LIBRARY_TEMPLATE_VERSION})")

set (LIBRARY_TEMPLATE_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (library_template-doxygen.cmake)

include (library_template-doxygen-layout.cmake)

enable_testing ()

if (LIBRARY_TEMPLATE_USER_DOC)
    add_subdirectory (doc_usr)
endif ()

if (LIBRARY_TEMPLATE_DEV_DOC)
    add_subdirectory (doc_dev)
endif ()

include (library_template-doxygen-package.cmake)
