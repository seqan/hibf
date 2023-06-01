cmake_minimum_required (VERSION 3.15)

if (NOT TARGET hibf_lib)
    add_library ("hibf_lib" STATIC hierarchical_interleaved_bloom_filter.cpp)

    target_link_libraries ("hibf_lib" PUBLIC "hibf::hibf" "hibf_layout")
endif ()

add_subdirectory (detail)
