cmake_minimum_required (VERSION 3.15)

if (NOT TARGET hibf_build)
    add_library ("hibf_build" STATIC
                 compute_kmers.cpp
                 construct_ibf.cpp
                 initialise_build_tree.cpp
                 insert_into_ibf.cpp
                 parse_chopper_pack_header.cpp
                 parse_chopper_pack_line.cpp
                 read_chopper_pack_file.cpp
                 update_content_node_data.cpp
                 update_header_node_data.cpp)

    target_link_libraries ("hibf_build" PUBLIC "hibf::hibf")
endif ()
