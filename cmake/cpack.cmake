# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
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
set (CPACK_SOURCE_IGNORE_FILES "\\\\.git($|/)")

include (CPack)
