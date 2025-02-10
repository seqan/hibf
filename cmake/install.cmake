# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set (HIBF_EXPORT_TARGETS "hibf")
if (TARGET cereal)
    list (APPEND HIBF_EXPORT_TARGETS "cereal")
endif ()
if (TARGET simde)
    list (APPEND HIBF_EXPORT_TARGETS "simde")
endif ()

# cmake-format: off
install (TARGETS ${HIBF_EXPORT_TARGETS}
         EXPORT hibf_targets
         INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         FRAMEWORK DESTINATION ${CMAKE_INSTALL_LIBDIR})
# cmake-format: on

install (DIRECTORY "${HIBF_HEADER_PATH}/hibf" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install (EXPORT hibf_targets
         NAMESPACE seqan::
         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hibf
         EXPORT_LINK_INTERFACE_LIBRARIES
         FILE hibf-targets.cmake)

include (CMakePackageConfigHelpers)
configure_package_config_file (cmake/hibf-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/hibf-config.cmake
                               INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hibf)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/hibf-config.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hibf)

set (version_file "${CMAKE_CURRENT_BINARY_DIR}/cmake/hibf-config-version.cmake")
write_basic_package_version_file (
    ${version_file}
    VERSION ${HIBF_VERSION}
    COMPATIBILITY SameMajorVersion)
install (FILES ${version_file} DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hibf)

install (FILES "${HIBF_SOURCE_DIR}/LICENSE.md" "${HIBF_SOURCE_DIR}/README.md"
         DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/hibf)
install (DIRECTORY "${HIBF_SOURCE_DIR}/LICENSES" DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/doc/hibf)
