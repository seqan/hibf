# ------------------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/hibf/blob/main/LICENSE.md
# ------------------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

# Finds all files relative to the `test_base_path_` which satisfy the given file pattern.
#
# Example:
# hibf_test_files (header_files "/hibf/include" "*.hpp;*.h")
#
# The variable `header_files` will contain:
#   hibf/alphabet/adaptation/all.hpp
#   hibf/alphabet/adaptation/char.hpp
#   hibf/alphabet/adaptation/concept.hpp
#   hibf/alphabet/adaptation/uint.hpp
#   hibf/alphabet/all.hpp
#   hibf/alphabet/dna5_detail.hpp
#   ....
macro (hibf_test_files VAR test_base_path_ extension_wildcards)
    # test_base_path is /home/.../hibf/test/
    get_filename_component (test_base_path "${test_base_path_}" ABSOLUTE)
    file (RELATIVE_PATH test_base_path_relative "${CMAKE_CURRENT_SOURCE_DIR}" "${test_base_path}")
    # ./ is a hack to deal with empty test_base_path_relative
    set (test_base_path_relative "./${test_base_path_relative}")
    # collect all cpp files
    set (${VAR} "")
    foreach (extension_wildcard ${extension_wildcards})
        file (GLOB_RECURSE test_files
              RELATIVE "${test_base_path}"
              "${test_base_path_relative}/${extension_wildcard}")
        list (APPEND ${VAR} ${test_files})
    endforeach ()

    unset (test_base_path)
    unset (test_base_path_relative)
endmacro ()
