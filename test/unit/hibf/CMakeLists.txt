# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

add_subdirectories ()

hibf_test (config_test.cpp)
hibf_test (hierarchical_interleaved_bloom_filter_test.cpp)
hibf_test (interleaved_bloom_filter_test.cpp)
hibf_test (print_test.cpp)
