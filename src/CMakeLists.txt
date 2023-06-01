cmake_minimum_required (VERSION 3.15)

if (NOT TARGET hibfh)
    add_library ("hibfh" STATIC hierarchical_interleaved_bloom_filter.cpp)

    target_link_libraries ("hibfh" PUBLIC "hibf::hibf")
endif ()
