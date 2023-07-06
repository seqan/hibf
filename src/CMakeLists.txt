set (HIBF_SOURCE_FILES
     hierarchical_interleaved_bloom_filter.cpp
     detail/layout/simple_binning.cpp
     detail/layout/execute.cpp
     detail/layout/output.cpp
     detail/layout/compute_fpr_correction.cpp
     detail/layout/compute_layout.cpp
     detail/layout/hierarchical_binning.cpp
     detail/sketch/toolbox.cpp
     detail/sketch/hyperloglog.cpp
     detail/build/initialise_build_tree.cpp
     detail/build/insert_into_ibf.cpp
     detail/build/compute_kmers.cpp
     detail/build/update_header_node_data.cpp
     detail/build/construct_ibf.cpp
     detail/build/update_content_node_data.cpp)

add_library (hibf STATIC ${HIBF_SOURCE_FILES})
add_library (seqan::hibf ALIAS hibf)
target_compile_definitions (hibf PUBLIC ${HIBF_DEFINITIONS})
target_compile_options (hibf PUBLIC ${HIBF_CXX_FLAGS})
target_link_options (hibf PUBLIC ${HIBF_CXX_FLAGS})
target_link_libraries (hibf PUBLIC "${HIBF_LIBRARIES}")

target_include_directories (
    hibf PUBLIC "$<BUILD_INTERFACE:${HIBF_HEADER_PATH}>" "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
                "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/hibf/submodules>")

target_include_directories (hibf SYSTEM PUBLIC "$<BUILD_INTERFACE:${HIBF_DEPENDENCY_HEADER_PATHS}>")
