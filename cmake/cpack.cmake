# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file describes how HIBF will be packaged.

cmake_minimum_required (VERSION 3.7)

set (CPACK_GENERATOR "TXZ")

set (CPACK_PACKAGE_VERSION "${HIBF_VERSION}")
set (CPACK_PACKAGE_VENDOR "seqan")
set (CPACK_PACKAGE_CHECKSUM "SHA256")
set (CPACK_PACKAGE_ICON "${HIBF_SOURCE_DIR}/test/documentation/seqan_logo.svg")
set (CPACK_RESOURCE_FILE_LICENSE "${HIBF_SOURCE_DIR}/LICENSE.md")
set (CPACK_RESOURCE_FILE_README "${HIBF_SOURCE_DIR}/README.md")

# Source Package
set (CPACK_SOURCE_GENERATOR "TXZ")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.git($|/)")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.github/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/\.vscode/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/build/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/submodules/")
list (APPEND CPACK_SOURCE_IGNORE_FILES "/cmake/CPM.cmake")

# Already being called on source package
if (NOT CPM_DOWNLOAD_LOCATION)
    set (CPM_DOWNLOAD_LOCATION "${HIBF_SOURCE_DIR}/cmake/CPM.cmake")
endif ()

configure_file ("${HIBF_SOURCE_DIR}/cmake/cpack_install.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/cpack_install.cmake"
                @ONLY)
set (CPACK_INSTALL_SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/cpack_install.cmake")

include (CPack)
