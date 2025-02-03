# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

add_subdirectories ()

hibf_test (bit_vector_test.cpp)
hibf_test (config_test.cpp)
hibf_test (counting_vector_test.cpp)
hibf_test (hierarchical_interleaved_bloom_filter_test.cpp)
hibf_test (insert_iterator_test.cpp)
hibf_test (interleaved_bloom_filter_test.cpp)
hibf_test (path_test.cpp)
hibf_test (print_test.cpp)
hibf_test (timer_test.cpp)

# If available, use Intel's Software Development Emulator (SDE) to emulate AVX-512 instructions.
# See https://www.intel.com/content/www/us/en/developer/articles/tool/software-development-emulator.html
find_program (SDE64_BIN sde64)

if (SDE64_BIN)
    hibf_test (counting_vector_avx512_test.cpp SDE "${SDE64_BIN}")
    target_compile_options (counting_vector_avx512_test PRIVATE -mavx512bw -mavx512f)

    hibf_test (interleaved_bloom_filter_avx512_test.cpp SDE "${SDE64_BIN}")
    target_compile_options (interleaved_bloom_filter_avx512_test PRIVATE -mavx512bw -mavx512f)
else ()
    hibf_test (counting_vector_avx512_test.cpp)
    hibf_test (interleaved_bloom_filter_avx512_test.cpp)
endif ()
