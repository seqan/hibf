# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# cereal
set (HIBF_CEREAL_VERSION 1.3.2)
CPMDeclarePackage (cereal
                   NAME cereal
                   VERSION ${HIBF_CEREAL_VERSION}
                   GITHUB_REPOSITORY USCiLab/cereal
                   SYSTEM TRUE
                   OPTIONS "JUST_INSTALL_CEREAL ON")
# simde
set (HIBF_SIMDE_VERSION 0.7.6)
CPMDeclarePackage (simde
                   NAME simde
                   VERSION ${HIBF_SIMDE_VERSION}
                   DOWNLOAD_ONLY YES
                   GITHUB_REPOSITORY simd-everywhere/simde)
# benchmark
set (HIBF_BENCHMARK_VERSION 1.8.2)
CPMDeclarePackage (benchmark
                   NAME benchmark
                   VERSION ${HIBF_BENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF")
# googletest
set (HIBF_GOOGLETEST_VERSION 1.14.0)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${HIBF_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF")
