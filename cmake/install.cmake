install (TARGETS hibf
         EXPORT hibf_targets
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         FRAMEWORK DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY "${HIBF_HEADER_PATH}/hibf" TYPE INCLUDE)
# install submodule header files, e.g. all external dependencies in /home/user/hibf/submodules/*,
# in /include/hibf/submodules/<submodule>/include
foreach (submodule_dir ${HIBF_DEPENDENCY_HEADER_PATHS})
    # e.g. submodule_dir: (1) /home/user/hibf/submodules/sdsl-lite/include or (2) /usr/include
    # strip /home/user/hibf/submodules/ and /include part.
    file (RELATIVE_PATH submodule "${HIBF_SOURCE_DIR}/submodules" "${submodule_dir}/..")
    # submodule is either a single module name, like sdsl-lite or a relative path to a folder ../../../usr
    # skip relative folders and only keep submodules that reside in the submodules folder
    if (NOT submodule MATCHES "^\\.\\.") # skip relative folders
        install (DIRECTORY "${submodule_dir}" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/hibf/submodules/${submodule}")
    endif ()
endforeach ()

install (EXPORT hibf_targets
         NAMESPACE seqan::
         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hibf
         FILE hibf-config.cmake)

include (CMakePackageConfigHelpers)
set (version_file "${CMAKE_CURRENT_BINARY_DIR}/cmake/hibf-config-version.cmake")
write_basic_package_version_file (
    ${version_file}
    VERSION ${HIBF_VERSION}
    COMPATIBILITY AnyNewerVersion)
install (FILES ${version_file} DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hibf)

install (FILES "${HIBF_SOURCE_DIR}/LICENSE.md" "${HIBF_SOURCE_DIR}/README.md" TYPE DOC)
