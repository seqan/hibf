# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.10...3.22)
project (hibf_test)

option (INSTALL_HIBF "" OFF)
enable_testing ()
add_subdirectory (../unit/ unit)
add_subdirectory (../performance/ performance)
add_subdirectory (../snippet/ snippet)
