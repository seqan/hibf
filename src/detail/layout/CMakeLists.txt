cmake_minimum_required (VERSION 3.15)

if (NOT TARGET hibf_layout)
    add_library ("hibf_layout" STATIC compute_fp_correction.cpp execute.cpp hierarchical_binning.cpp output.cpp
                                      simple_binning.cpp)

    target_link_libraries ("hibf_layout" PUBLIC "hibf::hibf")
endif ()
