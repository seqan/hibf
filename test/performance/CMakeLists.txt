# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10...3.22)
project (hibf_test_performance CXX)

include (../hibf-test.cmake)

set (HIBF_BENCHMARK_MIN_TIME
     "1"
     CACHE STRING "Set --benchmark_min_time= for each bechmark. Timings are unreliable in CI.")

macro (hibf_benchmark benchmark_cpp)
    file (RELATIVE_PATH benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${benchmark_cpp}")
    hibf_test_component (target "${benchmark}" TARGET_NAME)
    hibf_test_component (test_name "${benchmark}" TEST_NAME)

    add_executable (${target} ${benchmark_cpp})
    target_link_libraries (${target} hibf::test::performance)
    add_test (NAME "${test_name}" COMMAND ${target} "--benchmark_min_time=${HIBF_BENCHMARK_MIN_TIME}")

    unset (benchmark)
    unset (target)
    unset (test_name)
endmacro ()

hibf_require_benchmark ()

add_subdirectories ()
