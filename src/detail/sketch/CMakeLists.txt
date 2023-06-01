cmake_minimum_required (VERSION 3.15)

if (NOT TARGET hibf_sketch)
    add_library ("hibf_sketch" STATIC hyperloglog.cpp toolbox.cpp)

    target_link_libraries ("hibf_sketch" PUBLIC "hibf::hibf")
endif ()
