# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library_template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# Library_Template. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.4)

find_path (LIBRARY_TEMPLATE_MODULE_PATH "library_template-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system/")
list (APPEND CMAKE_MODULE_PATH "${LIBRARY_TEMPLATE_MODULE_PATH}")

include (library_template-config-version)

if (CMAKE_VERSION VERSION_LESS 3.12)
    project (library_template
             LANGUAGES CXX
             VERSION "${LIBRARY_TEMPLATE_PROJECT_VERSION}")
else ()
    project (library_template
             LANGUAGES CXX
             VERSION "${LIBRARY_TEMPLATE_PROJECT_VERSION}"
             DESCRIPTION "Library_Template -- the modern C++ library for sequence analysis" # since cmake 3.9
             HOMEPAGE_URL "https://github.com/seqan/library_template" # since cmake 3.12
    )
endif ()

find_package (Library_Template 3.0 REQUIRED HINTS ${LIBRARY_TEMPLATE_MODULE_PATH})

option (INSTALL_LIBRARY_TEMPLATE
        "Enable installation of library_template. (Projects embedding library_template may want to turn this OFF.)" ON)

if (INSTALL_LIBRARY_TEMPLATE)
    include (library_template-install)
    include (library_template-package)
endif ()
