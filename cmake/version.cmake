# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# -------------------------------------------------------------------------------------------------------------

file (STRINGS "${CMAKE_CURRENT_LIST_DIR}/../include/hibf/version.hpp" HIBF_VERSION_HPP
      REGEX "#define HIBF_VERSION_(MAJOR|MINOR|PATCH)")
string (REGEX REPLACE "#define HIBF_VERSION_(MAJOR|MINOR|PATCH) " "" HIBF_VERSION "${HIBF_VERSION_HPP}")
string (REGEX REPLACE ";" "." HIBF_VERSION "${HIBF_VERSION}")

file (STRINGS "${CMAKE_CURRENT_LIST_DIR}/../include/hibf/version.hpp" HIBF_RELEASE_CANDIDATE_HPP
      REGEX "#define HIBF_RELEASE_CANDIDATE ")
string (REGEX REPLACE "#define HIBF_RELEASE_CANDIDATE " "" HIBF_RELEASE_CANDIDATE_VERSION
                      "${HIBF_RELEASE_CANDIDATE_HPP}")
