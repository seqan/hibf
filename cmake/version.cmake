# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

file (STRINGS "${CMAKE_CURRENT_LIST_DIR}/../include/hibf/version.hpp" HIBF_VERSION_HPP
      REGEX "#define HIBF_VERSION_(MAJOR|MINOR|PATCH)")
string (REGEX REPLACE "#define HIBF_VERSION_(MAJOR|MINOR|PATCH) " "" HIBF_VERSION "${HIBF_VERSION_HPP}")
string (REGEX REPLACE ";" "." HIBF_VERSION "${HIBF_VERSION}")

file (STRINGS "${CMAKE_CURRENT_LIST_DIR}/../include/hibf/version.hpp" HIBF_RELEASE_CANDIDATE_HPP
      REGEX "#define HIBF_RELEASE_CANDIDATE ")
string (REGEX REPLACE "#define HIBF_RELEASE_CANDIDATE " "" HIBF_RELEASE_CANDIDATE_VERSION
                      "${HIBF_RELEASE_CANDIDATE_HPP}")
