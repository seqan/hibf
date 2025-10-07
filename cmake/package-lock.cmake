# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# cmake-format: off

# cereal
set (HIBF_CEREAL_VERSION 1.3.2 CACHE STRING "" FORCE)
CPMDeclarePackage (cereal
                   NAME cereal
                   VERSION ${HIBF_CEREAL_VERSION}
                   GITHUB_REPOSITORY USCiLab/cereal
                   SYSTEM TRUE
                   OPTIONS "JUST_INSTALL_CEREAL ON" "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# simde
set (HIBF_SIMDE_VERSION 0.8.2 CACHE STRING "" FORCE)
CPMDeclarePackage (simde
                   NAME simde
                   VERSION ${HIBF_SIMDE_VERSION}
                   GITHUB_REPOSITORY simd-everywhere/simde
                   DOWNLOAD_ONLY YES QUIET YES)
# benchmark
set (HIBF_BENCHMARK_VERSION 1.9.4 CACHE STRING "" FORCE)
CPMDeclarePackage (benchmark
                   NAME benchmark
                   VERSION ${HIBF_BENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# googletest
set (HIBF_GOOGLETEST_VERSION 1.17.0 CACHE STRING "" FORCE)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${HIBF_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING")
# doxygen-awesome
set (HIBF_DOXYGEN_AWESOME_VERSION 2.4.0 CACHE STRING "" FORCE)
CPMDeclarePackage (doxygen_awesome
                   NAME doxygen_awesome
                   VERSION ${HIBF_DOXYGEN_AWESOME_VERSION}
                   GITHUB_REPOSITORY jothepro/doxygen-awesome-css
                   DOWNLOAD_ONLY TRUE)

# cmake-format: on
