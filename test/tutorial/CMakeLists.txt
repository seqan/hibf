# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.16)
project (hibf_test_tutorial CXX)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

set (CPM_INDENT "CMake Package Manager CPM: ")
include (${CMAKE_CURRENT_LIST_DIR}/../../cmake/CPM.cmake)

# Tells CPM to use the version of hibf in the parent directory.
set (CPM_hibf_SOURCE "${CMAKE_CURRENT_LIST_DIR}/../..")
CPMAddPackage (NAME app_template
               GITHUB_REPOSITORY seqan/app-template
               GIT_TAG 718c88108bba2c04dc59855934feea9ac6f3a4a4
               OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING" "app-template_TEST OFF"
               EXCLUDE_FROM_ALL true)

add_executable ("hibf_tutorial" build.cpp count.cpp main.cpp search.cpp)
target_link_libraries ("hibf_tutorial" PRIVATE "app-template_interface")

FetchContent_Declare (
    download_data
    URL https://ftp.seqan.de/tutorial/hibf/tutorial_files.tar.gz
    URL_HASH SHA256=e07dd741196cffdc66ae4b0430ca64bdefd50c2af2aff652466585cbb56dceec)
FetchContent_Populate (download_data)

file (REMOVE "${download_data_SOURCE_DIR}/filenames.txt")
execute_process (COMMAND "${download_data_SOURCE_DIR}/get_filenames.sh"
                 WORKING_DIRECTORY "${download_data_SOURCE_DIR}"
                 OUTPUT_QUIET)

enable_testing ()
execute_process (COMMAND ${CMAKE_COMMAND} -E make_directory output)

add_test (NAME remove_old_output COMMAND ${CMAKE_COMMAND} -E remove output/index output/search output/count)
add_test (NAME build_index COMMAND hibf_tutorial build --input "${download_data_SOURCE_DIR}/filenames.txt" --output
                                   output/index --kmer 20)
add_test (NAME search_index COMMAND hibf_tutorial search --index output/index --query
                                    "${download_data_SOURCE_DIR}/reads.fastq" --output output/search --threshold 0.9)
add_test (NAME count_index COMMAND hibf_tutorial count --index output/index --query
                                   "${download_data_SOURCE_DIR}/reads.fastq" --output output/count --threshold 200)

set_tests_properties (remove_old_output build_index PROPERTIES FIXTURES_SETUP tutorial)
set_tests_properties (search_index count_index PROPERTIES FIXTURES_REQUIRED tutorial)
set_tests_properties (build_index PROPERTIES DEPENDS remove_old_output)
