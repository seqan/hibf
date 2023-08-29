set (HIBF_SOURCE_FILES
     hierarchical_interleaved_bloom_filter.cpp
     config.cpp
     layout/simple_binning.cpp
     layout/layout.cpp
     layout/execute.cpp
     layout/compute_fpr_correction.cpp
     layout/compute_layout.cpp
     layout/graph.cpp
     layout/hierarchical_binning.cpp
     sketch/toolbox.cpp
     sketch/hyperloglog.cpp
     build/insert_into_ibf.cpp
     build/compute_kmers.cpp
     build/construct_ibf.cpp)

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
