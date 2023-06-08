set (HIBF_SOURCE_FILES
     hierarchical_interleaved_bloom_filter.cpp
     detail/layout/simple_binning.cpp
     detail/layout/execute.cpp
     detail/layout/output.cpp
     detail/layout/compute_fp_correction.cpp
     detail/layout/hierarchical_binning.cpp
     detail/sketch/toolbox.cpp
     detail/sketch/hyperloglog.cpp
     detail/build/initialise_build_tree.cpp
     detail/build/insert_into_ibf.cpp
     detail/build/parse_chopper_pack_header.cpp
     detail/build/compute_kmers.cpp
     detail/build/read_chopper_pack_file.cpp
     detail/build/update_header_node_data.cpp
     detail/build/parse_chopper_pack_line.cpp
     detail/build/construct_ibf.cpp
     detail/build/update_content_node_data.cpp)

add_library (hibf STATIC ${HIBF_SOURCE_FILES})
target_compile_definitions (hibf PUBLIC ${HIBF_DEFINITIONS})
target_compile_options (hibf PUBLIC ${HIBF_CXX_FLAGS_LIST})
target_link_options (hibf PUBLIC ${HIBF_CXX_FLAGS_LIST})
target_link_libraries (hibf PUBLIC "${HIBF_LIBRARIES}")

target_include_directories (
    hibf PUBLIC "$<BUILD_INTERFACE:${HIBF_HEADER_PATH}>" "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
                "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hibf/submodules>")

target_include_directories (hibf SYSTEM PUBLIC "$<BUILD_INTERFACE:${HIBF_DEPENDENCY_INCLUDE_DIRS}>")

if (HIBF_HAS_LTO)
    target_compile_options (hibf PRIVATE ${HIBF_LTO_FLAGS})
endif ()

install (TARGETS hibf
         EXPORT hibf_targets
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
         FRAMEWORK DESTINATION ${CMAKE_INSTALL_LIBDIR})

install (DIRECTORY "${HIBF_HEADER_PATH}/hibf" TYPE INCLUDE)
# install submodule header files, e.g. all external dependencies in /home/user/hibf/submodules/*,
# in /include/hibf/submodules/<submodule>/include
foreach (submodule_dir ${HIBF_DEPENDENCY_INCLUDE_DIRS})
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
