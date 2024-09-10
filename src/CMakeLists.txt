# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set (HIBF_SOURCE_FILES
     hierarchical_interleaved_bloom_filter.cpp
     interleaved_bloom_filter.cpp
     config.cpp
     layout/simple_binning.cpp
     layout/layout.cpp
     layout/compute_fpr_correction.cpp
     layout/compute_layout.cpp
     layout/compute_relaxed_fpr_correction.cpp
     sketch/compute_sketches.cpp
     layout/graph.cpp
     layout/hierarchical_binning.cpp
     misc/print.cpp
     sketch/toolbox.cpp
     sketch/hyperloglog.cpp
     sketch/minhashes.cpp
     build/insert_into_ibf.cpp
     build/compute_kmers.cpp
     build/construct_ibf.cpp)

if (NOT TARGET hibf)
    add_library (hibf STATIC ${HIBF_SOURCE_FILES})
    add_library (seqan::hibf ALIAS hibf)
    target_compile_definitions (hibf PUBLIC ${HIBF_DEFINITIONS})
    target_compile_options (hibf PUBLIC ${HIBF_CXX_FLAGS})
    target_link_options (hibf PUBLIC ${HIBF_CXX_FLAGS})
    target_link_libraries (hibf PUBLIC "${HIBF_LIBRARIES}")
    target_link_libraries (hibf PUBLIC cereal::cereal simde::simde OpenMP::OpenMP_CXX)

    target_include_directories (hibf PUBLIC "$<BUILD_INTERFACE:${HIBF_HEADER_PATH}>"
                                            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
endif ()
