# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# This file describes how Library_Template will be packaged.

cmake_minimum_required (VERSION 3.7)

set (CPACK_GENERATOR "TXZ")

set (CPACK_PACKAGE_VERSION "${LIBRARY_TEMPLATE_VERSION}")
set (CPACK_PACKAGE_VENDOR "seqan")
set (CPACK_PACKAGE_CHECKSUM "SHA256")
set (CPACK_PACKAGE_ICON "${LIBRARY_TEMPLATE_CLONE_DIR}/test/documentation/seqan_logo.svg")
set (CPACK_RESOURCE_FILE_LICENSE "${LIBRARY_TEMPLATE_CLONE_DIR}/LICENSE.md")
set (CPACK_RESOURCE_FILE_README "${LIBRARY_TEMPLATE_CLONE_DIR}/README.md")

# Source Package
set (CPACK_SOURCE_GENERATOR "TXZ")
set (CPACK_SOURCE_IGNORE_FILES "\\\\.git($|/)")

include (CPack)
