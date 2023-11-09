# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# 1: Check via pkg-config. Debian package provides a pc file.
# 2: Via CPM

macro (hibf_define_simde SIMDE_INCLUDE_DIRECTORY)
    add_library (simde INTERFACE)
    add_library (simde::simde ALIAS simde)
    target_include_directories (simde INTERFACE "$<BUILD_INTERFACE:${SIMDE_INCLUDE_DIRECTORY}>")
    set_target_properties (simde PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                                            $<TARGET_PROPERTY:simde,INTERFACE_INCLUDE_DIRECTORIES>)
endmacro ()

if (CPM_USE_LOCAL_PACKAGES)
    find_package (PkgConfig QUIET)
    if (PKG_CONFIG_FOUND)
        pkg_check_modules (simde QUIET simde>=${HIBF_SIMDE_VERSION})
    endif ()
endif ()

if (simde_FOUND)
    hibf_define_simde ("${simde_INCLUDEDIR}")
    cpm_message (STATUS "${CPM_INDENT} Using local package simde@${simde_VERSION}")
else ()
    CPMGetPackage (simde)

    if (simde_ADDED)
        hibf_define_simde ("${simde_SOURCE_DIR}")
    endif ()
endif ()
