# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/library-template/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10...3.22)
project (library_template_test_performance CXX)

include (../library_template-test.cmake)

set (LIBRARY_TEMPLATE_BENCHMARK_MIN_TIME
     "1"
     CACHE STRING "Set --benchmark_min_time= for each bechmark. Timings are unreliable in CI.")

macro (library_template_benchmark benchmark_cpp)
    file (RELATIVE_PATH benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${benchmark_cpp}")
    library_template_test_component (target "${benchmark}" TARGET_NAME)
    library_template_test_component (test_name "${benchmark}" TEST_NAME)

    add_executable (${target} ${benchmark_cpp})
    target_link_libraries (${target} library_template::test::performance)
    add_test (NAME "${test_name}" COMMAND ${target} "--benchmark_min_time=${LIBRARY_TEMPLATE_BENCHMARK_MIN_TIME}")

    unset (benchmark)
    unset (target)
    unset (test_name)
endmacro ()

library_template_require_ccache ()
library_template_require_benchmark ()

add_subdirectories ()
